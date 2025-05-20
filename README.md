# Mapperia Map Generator

Generate crisp, laser-ready SVG and preview PNG maps from OpenStreetMap data. Choose layers (buildings, streets, greens, water, railways) and export for preview or laser cutting.

## Tech Stack

- **Backend:** FastAPI, svgwrite, cairosvg, overpy, geopandas, shapely, pyproj
- **Frontend:** Next.js (React, TypeScript, Tailwind)

## Features

- Exports: **SVG (laser-ready, thin outlines), PNG (adaptive resolution)**
- Adaptive scaling for all map sizes

## Usage

### Backend

```bash
cd backend
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
uvicorn src.api:app --reload
```

### Frontend

```bash
cd frontend
npm install
npm run dev
```
