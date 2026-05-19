// One-page vMPT summary slide. Run with:
//   NODE_PATH=$(npm root -g) node build_vmpt_slide.js
const pptxgen = require("pptxgenjs");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const sharp = require("sharp");
const {
  FaTelescope, FaCrosshairs, FaMousePointer, FaFileExport,
  FaCheckCircle, FaCode, FaSatellite,
} = require("react-icons/fa");
const { MdGridOn, MdHighlight, MdRouter } = require("react-icons/md");
const { HiOutlineSparkles } = require("react-icons/hi");

function renderIconSvg(IconComponent, color, size) {
  return ReactDOMServer.renderToStaticMarkup(
    React.createElement(IconComponent, { color, size: String(size) })
  );
}
async function iconPng(IconComponent, color, size = 256) {
  const svg = renderIconSvg(IconComponent, color, size);
  const buf = await sharp(Buffer.from(svg)).png().toBuffer();
  return "image/png;base64," + buf.toString("base64");
}

// ---- Color palette: warm cream paper · authentic-dark instrument viewport ----
const C = {
  bg:       "F7F4ED",   // warm cream (slide background)
  panel:    "FFFFFF",   // white cards / header
  panelAlt: "F0EDE3",   // alt panel (subtle warm tint)
  border:   "D8D1BF",   // soft sand border
  text:     "1A2C4E",   // deep navy (primary text)
  muted:    "5C6B82",   // muted slate (captions)
  accent:   "C58A00",   // deep amber (titles / stats)
  accent2:  "2E5BBF",   // MSA blue (deepened for light bg)
  red:      "D63D3D",   // open-shutter red (slightly deeper)
  green:    "2E9B3F",   // pointing green (legible on light bg)
  orange:   "E26A00",   // spec overlap (deeper orange)
  // Mock-viewport-only (kept dark — sky is authentically dark)
  vpBg:     "0A1428",
  vpStar:   "FFFFFF",
  vpGalaxy: "FFD080",
  vpBlue:   "3A7EFF",
  vpRed:    "FF5454",
  vpGreen:  "7FFF4A",
  vpOrange: "FF8A1F",
};

(async () => {
  const pres = new pptxgen();
  pres.layout = "LAYOUT_WIDE"; // 13.3" x 7.5"
  pres.author = "Sunfeng Wu";
  pres.title  = "vMPT — visual MSA Planning Tool";

  const s = pres.addSlide();
  s.background = { color: C.bg };

  // Deterministic RNG (no slide-wide starfield on the light bg — keep it clean).
  const rand = (() => { let x = 12345; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff, (x >>> 0) / 0x7fffffff); })();

  // ─── Header bar ─────────────────────────────────────────────────────────────
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: 0, w: 13.3, h: 1.05,
    fill: { color: C.panel }, line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {  // thin separator below header
    x: 0, y: 1.05, w: 13.3, h: 0.03,
    fill: { color: C.accent }, line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {  // accent stripe along left edge of header
    x: 0, y: 0, w: 0.12, h: 1.05,
    fill: { color: C.accent }, line: { type: "none" },
  });
  // Logo-style favicon mark (mini MSA grid)
  const logoX = 0.45, logoY = 0.18, logoW = 0.7;
  s.addShape(pres.shapes.RECTANGLE, {
    x: logoX, y: logoY, w: logoW, h: logoW,
    fill: { color: C.vpBg }, line: { color: C.accent2, width: 1 },
  });
  // 4 quadrants of MSA inside the logo
  const qPad = 0.07, qW = (logoW - qPad * 3) / 2;
  for (let r = 0; r < 2; r++) for (let c = 0; c < 2; c++) {
    s.addShape(pres.shapes.RECTANGLE, {
      x: logoX + qPad + c * (qW + qPad), y: logoY + qPad + r * (qW + qPad),
      w: qW, h: qW,
      fill: { type: "none" }, line: { color: C.accent2, width: 0.75 },
    });
  }
  // lime "pointing" cross
  s.addShape(pres.shapes.LINE, {
    x: logoX + logoW/2 - 0.10, y: logoY + logoW/2,
    w: 0.20, h: 0, line: { color: C.vpGreen, width: 2 },
  });
  s.addShape(pres.shapes.LINE, {
    x: logoX + logoW/2, y: logoY + logoW/2 - 0.10,
    w: 0, h: 0.20, line: { color: C.vpGreen, width: 2 },
  });

  s.addText("vMPT", {
    x: 1.35, y: 0.05, w: 4.0, h: 0.60, margin: 0,
    fontFace: "Georgia", fontSize: 44, bold: true, color: C.text, charSpacing: 2,
  });
  s.addText("visual MSA Planning Tool", {
    x: 1.35, y: 0.62, w: 5.5, h: 0.40, margin: 0,
    fontFace: "Georgia", fontSize: 17, italic: true, color: C.muted,
  });
  s.addText("JWST / NIRSpec  •  interactive shutter picking + APT export", {
    x: 6.0, y: 0.36, w: 7.1, h: 0.45, margin: 0,
    fontFace: "Calibri", fontSize: 15, color: C.muted, align: "right",
  });

  // ─── Mock app viewport (left half) ─────────────────────────────────────────
  // The sky is intentionally dark — this is the actual app's image canvas.
  const vpX = 0.45, vpY = 1.30, vpW = 5.8, vpH = 5.35;
  s.addShape(pres.shapes.RECTANGLE, {
    x: vpX, y: vpY, w: vpW, h: vpH,
    fill: { color: C.vpBg }, line: { color: C.border, width: 1 },
  });
  // viewport title bar
  s.addShape(pres.shapes.RECTANGLE, {
    x: vpX, y: vpY, w: vpW, h: 0.36,
    fill: { color: "152846" }, line: { type: "none" },
  });
  s.addText("Image · MSA overlay · target picks", {
    x: vpX + 0.15, y: vpY + 0.03, w: vpW - 0.3, h: 0.3, margin: 0,
    fontFace: "Calibri", fontSize: 12, color: "C0CCE0", italic: true,
  });

  // Sky background (subtle starfield inside the viewport)
  for (let i = 0; i < 70; i++) {
    const sz = 0.02 + rand() * 0.05;
    s.addShape(pres.shapes.OVAL, {
      x: vpX + 0.1 + rand() * (vpW - 0.2),
      y: vpY + 0.5 + rand() * (vpH - 0.7),
      w: sz, h: sz,
      fill: { color: C.vpStar, transparency: 50 + Math.floor(rand() * 45) },
      line: { type: "none" },
    });
  }
  // a few brighter "galaxies"
  for (let i = 0; i < 12; i++) {
    const sz = 0.10 + rand() * 0.10;
    s.addShape(pres.shapes.OVAL, {
      x: vpX + 0.3 + rand() * (vpW - 0.6),
      y: vpY + 0.7 + rand() * (vpH - 1.0),
      w: sz, h: sz,
      fill: { color: C.vpGalaxy, transparency: 60 },
      line: { type: "none" },
    });
  }

  // MSA: 4 quadrants rotated ~-20° (illusion via slight shift)
  // We approximate by drawing 4 rectangles in a 2x2 with small gap.
  const msaCx = vpX + vpW/2, msaCy = vpY + 0.45 + (vpH - 0.45)/2 + 0.1;
  const qSide = 1.55, qGap = 0.05;
  const quadPositions = [
    { x: msaCx - qSide - qGap/2, y: msaCy - qSide - qGap/2 },
    { x: msaCx + qGap/2,         y: msaCy - qSide - qGap/2 },
    { x: msaCx - qSide - qGap/2, y: msaCy + qGap/2 },
    { x: msaCx + qGap/2,         y: msaCy + qGap/2 },
  ];
  // First, draw faint orange "spec overlap" bands across some rows
  const overlapRows = [
    { qIdx: 0, rowFrac: 0.30 },
    { qIdx: 1, rowFrac: 0.30 },
    { qIdx: 2, rowFrac: 0.55 },
    { qIdx: 3, rowFrac: 0.55 },
    { qIdx: 0, rowFrac: 0.70 },
    { qIdx: 1, rowFrac: 0.70 },
  ];
  // Band spans BOTH quadrants in the same row → use full MSA width
  const seenRows = new Set();
  for (const ov of overlapRows) {
    const key = `${Math.floor(ov.qIdx / 2)}-${ov.rowFrac}`;
    if (seenRows.has(key)) continue;
    seenRows.add(key);
    const top = (ov.qIdx < 2) ? quadPositions[0] : quadPositions[2];
    const bandY = top.y + ov.rowFrac * qSide - 0.06;
    s.addShape(pres.shapes.RECTANGLE, {
      x: msaCx - qSide - qGap/2, y: bandY,
      w: qSide * 2 + qGap, h: 0.12,
      fill: { color: C.vpOrange, transparency: 75 }, line: { type: "none" },
    });
  }
  // MSA quadrant outlines
  for (const q of quadPositions) {
    s.addShape(pres.shapes.RECTANGLE, {
      x: q.x, y: q.y, w: qSide, h: qSide,
      fill: { type: "none" },
      line: { color: C.vpBlue, width: 1.5 },
    });
  }
  // Operable-shutter dots (silver edge stipple) — represent with a faint dot grid
  const dotsPerSide = 9;
  for (const q of quadPositions) {
    for (let i = 0; i < dotsPerSide; i++) {
      for (let j = 0; j < dotsPerSide; j++) {
        const x = q.x + 0.10 + i * ((qSide - 0.20) / (dotsPerSide - 1));
        const y = q.y + 0.10 + j * ((qSide - 0.20) / (dotsPerSide - 1));
        s.addShape(pres.shapes.OVAL, {
          x: x - 0.015, y: y - 0.015, w: 0.03, h: 0.03,
          fill: { color: "C0C0C0", transparency: 75 }, line: { type: "none" },
        });
      }
    }
  }
  // Open shutters (red 3-shutter slitlets)
  const slitlets = [
    { q: 0, ix: 2, iy: 4 }, { q: 0, ix: 5, iy: 2 },
    { q: 1, ix: 3, iy: 5 }, { q: 1, ix: 6, iy: 3 },
    { q: 2, ix: 1, iy: 6 }, { q: 2, ix: 4, iy: 2 },
    { q: 3, ix: 2, iy: 4 }, { q: 3, ix: 6, iy: 6 }, { q: 3, ix: 5, iy: 1 },
  ];
  for (const sl of slitlets) {
    const q = quadPositions[sl.q];
    const cellW = (qSide - 0.20) / (dotsPerSide - 1);
    const cx = q.x + 0.10 + sl.ix * cellW;
    const cy = q.y + 0.10 + sl.iy * cellW;
    // 3 stacked shutters
    for (let k = -1; k <= 1; k++) {
      s.addShape(pres.shapes.RECTANGLE, {
        x: cx - 0.05, y: cy + k * 0.08 - 0.035, w: 0.10, h: 0.07,
        fill: { color: C.vpRed, transparency: 50 }, line: { color: C.vpRed, width: 0.5 },
      });
    }
  }
  // Stuck-open shutters (thick dark-red outline, hollow)
  const stuck = [{ q: 0, ix: 8, iy: 0 }, { q: 3, ix: 0, iy: 8 }];
  for (const sk of stuck) {
    const q = quadPositions[sk.q];
    const cellW = (qSide - 0.20) / (dotsPerSide - 1);
    const cx = q.x + 0.10 + sk.ix * cellW;
    const cy = q.y + 0.10 + sk.iy * cellW;
    s.addShape(pres.shapes.RECTANGLE, {
      x: cx - 0.06, y: cy - 0.04, w: 0.12, h: 0.08,
      fill: { type: "none" }, line: { color: "B30000", width: 2 },
    });
  }
  // Pointing cross at center
  s.addShape(pres.shapes.LINE, {
    x: msaCx - 0.18, y: msaCy, w: 0.36, h: 0,
    line: { color: C.vpGreen, width: 2.5 },
  });
  s.addShape(pres.shapes.LINE, {
    x: msaCx, y: msaCy - 0.18, w: 0, h: 0.36,
    line: { color: C.vpGreen, width: 2.5 },
  });

  // Mini legend along the bottom of the viewport
  const legY = vpY + vpH - 0.35;
  const legendItems = [
    { c: C.vpBlue,   label: "MSA quadrants" },
    { c: C.vpRed,    label: "Open shutters" },
    { c: C.vpOrange, label: "Spec overlap" },
    { c: C.vpGreen,  label: "Pointing"      },
  ];
  let legX = vpX + 0.15;
  for (const it of legendItems) {
    s.addShape(pres.shapes.RECTANGLE, {
      x: legX, y: legY + 0.06, w: 0.12, h: 0.12,
      fill: { color: it.c }, line: { type: "none" },
    });
    s.addText(it.label, {
      x: legX + 0.20, y: legY, w: 1.30, h: 0.24, margin: 0,
      fontFace: "Calibri", fontSize: 11, color: "C0CCE0", valign: "middle",
    });
    legX += 1.40;
  }

  // ─── Right column: feature cards + workflow ────────────────────────────────
  const rX = 6.55, rW = 6.40;
  const features = [
    {
      icon: FaSatellite, color: C.accent2,
      title: "Real instrument geometry",
      body: "Per-shutter (q,s,d) V2/V3 from pysiaf · MSA-on-sky for any V3 PA · live operability + stuck-open map from CRDS msaoper.",
    },
    {
      icon: FaMousePointer, color: C.green,
      title: "Hand-pick · snap · undo",
      body: "Click → snap to nearest operable shutter; shift-click moves pointing; auto 3-shutter slitlets on targets; cyan double-tap highlights.",
    },
    {
      icon: MdGridOn, color: C.orange,
      title: "Spec-overlap aware",
      body: "Per-grating V2 half-extent (PRISM 35″ → G*H 500″). Cross-quadrant Q1↔Q3 / Q2↔Q4 detector pairing. Stuck-open contributes overlap.",
    },
    {
      icon: FaFileExport, color: C.accent,
      title: "APT-ready export",
      body: "Pure MPT plan JSON (matches APT reference schema field-for-field) + eMPT shutter_mask.csv + observed_targets.cat + pointing_summary.txt.",
    },
  ];
  const cardX = rX, cardW = rW, cardH = 1.12, cardGapY = 0.14;
  const cardY0 = 1.30;

  for (let i = 0; i < features.length; i++) {
    const f = features[i];
    const y = cardY0 + i * (cardH + cardGapY);
    // Card background
    s.addShape(pres.shapes.RECTANGLE, {
      x: cardX, y, w: cardW, h: cardH,
      fill: { color: C.panel }, line: { color: C.border, width: 1 },
    });
    // Left accent bar
    s.addShape(pres.shapes.RECTANGLE, {
      x: cardX, y, w: 0.08, h: cardH,
      fill: { color: f.color }, line: { type: "none" },
    });
    // Icon in a small tinted circle (pale fill of the feature color)
    const iconData = await iconPng(f.icon, "#" + f.color, 256);
    const iconBox = 0.52;
    s.addShape(pres.shapes.OVAL, {
      x: cardX + 0.22, y: y + (cardH - iconBox - 0.22) / 2,
      w: iconBox + 0.22, h: iconBox + 0.22,
      fill: { color: f.color, transparency: 85 }, line: { color: f.color, width: 1.5 },
    });
    s.addImage({
      data: iconData,
      x: cardX + 0.33, y: y + (cardH - iconBox) / 2,
      w: iconBox, h: iconBox,
    });
    // Title + body
    s.addText(f.title, {
      x: cardX + 1.15, y: y + 0.13, w: cardW - 1.25, h: 0.38, margin: 0,
      fontFace: "Georgia", fontSize: 18, bold: true, color: C.text,
    });
    s.addText(f.body, {
      x: cardX + 1.15, y: y + 0.52, w: cardW - 1.25, h: cardH - 0.58, margin: 0,
      fontFace: "Calibri", fontSize: 13, color: C.muted, paraSpaceAfter: 0,
    });
  }

  // ─── Footer: tech stack + stats ────────────────────────────────────────────
  const footY = 6.78;
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: footY, w: 13.3, h: 0.72,
    fill: { color: C.panelAlt }, line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: footY, w: 13.3, h: 0.04,
    fill: { color: C.accent }, line: { type: "none" },
  });

  // Stats cluster on the left
  const stats = [
    { num: "249,660", lbl: "shutters / image" },
    { num: "60+",     lbl: "tests passing"     },
    { num: "<70 ms",  lbl: "redraw at zoom"    },
    { num: "MIT",     lbl: "license · open"    },
  ];
  let stX = 0.45;
  for (const st of stats) {
    s.addText(st.num, {
      x: stX, y: footY + 0.06, w: 1.7, h: 0.34, margin: 0,
      fontFace: "Georgia", fontSize: 20, bold: true, color: C.accent,
    });
    s.addText(st.lbl, {
      x: stX, y: footY + 0.38, w: 1.7, h: 0.24, margin: 0,
      fontFace: "Calibri", fontSize: 11, color: C.muted,
    });
    stX += 1.7;
  }

  // Stack chip on the right
  s.addText(
    "Bokeh · pysiaf · astropy · jwst_gtvt · CRDS msaoper",
    {
      x: 6.9, y: footY + 0.19, w: 6.15, h: 0.32, margin: 0,
      fontFace: "Consolas", fontSize: 12, color: C.muted, align: "right",
    }
  );

  await pres.writeFile({ fileName: "/Users/sunfengwu/nirspec/vmpt_summary.pptx" });
  console.log("wrote vmpt_summary.pptx");
})();
