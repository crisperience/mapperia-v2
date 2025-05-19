# Mapperia Map Generator

This tool generates SVG and PNG maps based on OpenStreetMap data. You can select layers (buildings, streets, greens, water, railways) and export the map for preview or laser cutting.

## Getting Started

### Backend

```bash
cd backend
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
uvicorn src.api:app --reload
```

The backend will be available at http://localhost:8000

### Frontend

```bash
cd frontend
npm install
npm run dev
```

The frontend will be available at http://localhost:3000
