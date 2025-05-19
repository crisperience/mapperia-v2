# Mapperia Map Generator

## What Does It Do?

- **Generates custom maps** for laser engraving or printing.
- **Web interface** lets you select areas, styles, and options.
- **Backend** processes your choices and outputs ready-to-use files.
- **Supports custom data layers** (buildings, streets, parks, water, etc.).

---

## How It Works?

1. **You open the web app** in your browser.
2. **Choose a location and map style.**
3. **Click generate.**
4. **The backend creates files** (SVG / PNG) for you to download and use.

---

## Quickstart for Developers

### Prerequisites

- Python 3.10+
- Node.js 18+
- (Recommended) Use a virtual environment for Python

### 1. Backend (Python)

```
cd backend
pip install -r requirements.txt
uvicorn src.api:app --reload
```

### 2. Frontend (Next.js)

```
cd frontend
npm install
npm run dev
```

- The frontend will run on [http://localhost:3000](http://localhost:3000)
- The backend will run on [http://localhost:8000](http://localhost:8000)

---

## Project Structure

- `backend/` — Python backend (map generation, data processing)
  - `src/` — Main backend code
  - `gpkg/` — Map data files (GeoPackage: buildings, streets, parks, water, etc.)
  - `output/` — Generated map files
- `frontend/` — Next.js web app (user interface)

---
