# 7.1-Sequence-Alignment-App

A Flask-based sequence alignment app for comparing DNA or protein sequences using global and local alignment algorithms from the package.

## Features
- Compare two sequences side by side with global or local alignment
- Choose between DNA and protein workflows
- Use guided scoring profiles for simple, balanced, strict, permissive, and protein-based comparisons
- Load example presets for SARS-CoV and SARS-CoV-2 spike sequences fetched live from NCBI
- View alignment metrics including matches, mismatches, gaps, and identity

## Project structure
- app.py: Flask web app entry point
- align_algorithms.py: global and local alignment implementations
- metrics.py: alignment summary metrics
- diagnostics.py: explanatory labels and interpretation text
- scoring.py: scoring and conservative-substitution logic
- templates/: HTML templates for the web interface
- static/: CSS styling for the app
- tests/: regression tests for alignment metrics

## Run locally
1. Create and activate a virtual environment if needed.
2. Install Flask if it is not already available:
   ```bash
   pip install Flask
   ```
3. Start the application:
   ```bash
   python app.py
   ```
4. Open http://127.0.0.1:5000/ in a browser.

## Notes
- DNA comparisons now count all mismatches correctly, even for short sequences such as CTCTCT vs CTCAGA.
- Protein comparisons still support conservative-substitution tracking where appropriate.
