# Molecule viewer

A browser-based viewer for molecules and LAMMPS simulation structures, with custom
colors and labels, a Studio renderer (materials, lights, shadows, ambient
occlusion), and high-resolution image, GIF, and MP4 exports.

Live at **https://sahmed73.github.io/molecule-viewer/**.

## How it is published

GitHub Pages serves this repository from the root of the `main` branch. Because the
account's user site is `sahmed73.github.io`, a repository named `molecule-viewer` is
published at `sahmed73.github.io/molecule-viewer/`, inside the website. Every push to
`main` goes live within a minute or two; nothing has to change in the website
repository. `.nojekyll` makes Pages serve the files as they are.

The page shares the website's look without copying it:

- `index.html` loads the website's stylesheet from `/assets/main.css`.
- `assets/site-chrome.mjs` replaces the page's simple header and footer with the
  website's current ones, read from `/assets/site-chrome.json`. The website builds
  that file from its own templates, so navigation changes there appear here too.
- `assets/molecule-viewer.css` holds the viewer's own styles, plus fallback base styles
  used only when the website's stylesheet is not available.

## Run locally

Serve the repository root over HTTP (ES modules do not load from `file://`):

```sh
python3 -m http.server 8000
```

Open `http://localhost:8000/`. Locally the page uses the fallback header and styles,
and the SMILES box can only use PubChem, because the RDKit service accepts requests
from the website only.

## Tests

```sh
node scripts/test_molecule_viewer.mjs
```

The checks cover the LAMMPS parser, styles, export sizes, video bitrates, Studio
geometry, GIF encoding, and that every example in the dropdown has a valid MOL file.

## How it works

- Structure files (LAMMPS `.data`, MOL/SDF, PDB, XYZ) are read in the browser and never
  uploaded.
- The Examples dropdown lists the MOL files in `assets/examples/`; each `<option>` in
  `index.html` names its file. Lubricant examples: PAO, ZDDP (dibutyl, with an
  illustrative Zn–S arrangement), ZDDP (sec-hexyl), and ADPA. Well-known molecules
  (water, caffeine, aspirin, glucose, ATP, cholesterol, C60, and more) use PubChem
  SMILES with RDKit 3D coordinates. `scripts/make_examples.py` rebuilds every example
  except PAO, ZDDP (dibutyl), and ADPA, and checks each one's stereochemistry against
  its PubChem InChIKey (needs `pip install rdkit`).
- SMILES requests a 3D structure from PubChem, then falls back to an RDKit service on
  Render. That service's code and `render.yaml` live in the website repository
  (`rdkit-api/`), and its URL is set in `data-rdkit-api` in `index.html`.
- [3Dmol.js](https://3dmol.csb.pitt.edu/) 2.0.4 (from cdnjs) draws the standard view
  and handles loading, picking, and rotation.
- The optional Studio renderer (`assets/molecule-studio.mjs`) redraws the same view
  with [three.js](https://threejs.org/) r170 from jsDelivr (see the import map in
  `index.html`). `assets/molecule-scene.mjs` builds its geometry the same way 3Dmol
  draws each style, and `assets/molecule-labels.mjs` draws the labels for both
  renderers.
- Exports render the 3D canvas itself at the chosen size, and labels scale with it.
  Rotation videos are encoded frame by frame with WebCodecs and packaged as MP4 by
  the bundled `mp4-muxer` 5.2.2; browsers without WebCodecs fall back to real-time
  recording. GIFs use the bundled `gifenc`.

## License

MIT (see `LICENSE`): anyone may use, change, and share the code, including in other
projects, as long as they keep the copyright notice. The bundled third-party files
below keep their own licenses.

## Third-party code

- `assets/mp4-muxer.mjs`: mp4-muxer 5.2.2, MIT license (`assets/mp4-muxer.LICENSE.md`).
- `assets/gifenc.mjs`: gifenc, MIT license (`assets/gifenc.LICENSE.md`).
- 3Dmol.js and three.js load from public CDNs at pinned versions.

## History

The viewer started inside the website repository
([sahmed73/sahmed73.github.io](https://github.com/sahmed73/sahmed73.github.io)). Its
commits were carried over here when it moved to its own repository.
