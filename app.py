import os
from pathlib import Path
import json
import uvicorn
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
import config.base_models as base_models
from manager import Manager
from utils.presets import get_preset_map, load_presets


BASE_DIR = Path(__file__).resolve().parent
print(str(Path(BASE_DIR, "templates")))
TEMPLATES = Jinja2Templates(directory=str(Path(BASE_DIR, "templates")))


app = FastAPI()
app.mount(
    "/static",
    StaticFiles(directory=str(Path(BASE_DIR, "static"))),
    name="static",
)


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    presets = [{key: value for key, value in preset.items() if key != "path"} for preset in load_presets(BASE_DIR)]
    context = {"request": request, "presets": presets}
    return TEMPLATES.TemplateResponse("index.html", context)


@app.post("/find-motif", response_model=base_models.FindMotifResponse)
async def find_motif(request_data: base_models.FindMotif):
    manager = Manager(request_data)
    response = manager.find_motif()
    return response


@app.get("/preset/{preset_slug}")
async def get_preset(preset_slug: str):
    preset = get_preset_map(BASE_DIR).get(preset_slug)
    if not preset:
        raise HTTPException(status_code=404, detail="Preset not found")
    preset_data = json.loads(Path(preset["path"]).read_text())
    if "DNA" not in preset_data:
        raise HTTPException(status_code=422, detail="Preset is missing DNA data")
    return preset_data


if __name__ == "__main__":
    port = int(os.getenv("PORT", 8080))
    uvicorn.run(app, port=port, log_level="info", host="0.0.0.0", access_log=False)
