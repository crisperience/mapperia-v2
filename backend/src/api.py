import os
import traceback

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel

from .main import generate_map

app = FastAPI()


class GenerateMapRequest(BaseModel):
    lat: float
    lon: float
    width: float
    height: float
    formats: list[str]
    style: str
    location: str
    layers: list[str]


@app.post("/api/generate-map")
async def generate_map_endpoint(payload: GenerateMapRequest):
    """Generate a map based on the provided parameters.

    Args:
        payload (GenerateMapRequest): Map generation parameters.

    Returns:
        JSONResponse: URLs to generated files or error message.
    """
    try:
        result = generate_map(
            lat=payload.lat,
            lon=payload.lon,
            width=payload.width,
            height=payload.height,
            formats=payload.formats,
            style=payload.style,
            location=payload.location,
            layers=payload.layers,
        )
        return JSONResponse(content=result)
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"error": str(e), "traceback": traceback.format_exc()},
        )


@app.get("/api/output/{filename}")
def get_output_file(filename: str):
    """Serve generated map files for download.

    Args:
        filename (str): Name of the file to serve.

    Returns:
        FileResponse: The requested file or 404 if not found.
    """
    output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output"))
    file_path = os.path.join(output_dir, filename)
    if not os.path.isfile(file_path):
        raise HTTPException(status_code=404, detail="File not found")
    # Force download for PNG and SVG
    if filename.lower().endswith((".png", ".svg")):
        return FileResponse(
            file_path, media_type="application/octet-stream", filename=filename
        )
    return FileResponse(file_path)
