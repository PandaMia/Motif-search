import os
from pathlib import Path
import uvicorn
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
import config.base_models as base_models
from manager import Manager


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
    context = {"request": request}
    return TEMPLATES.TemplateResponse("index.html", context)


@app.post("/find-motif", response_model=base_models.FindMotifResponse)
async def find_motif(request_data: base_models.FindMotif):
    manager = Manager(request_data)
    response = manager.find_motif()
    return response


@app.get("/preset/dosr", response_class=PlainTextResponse)
async def preset_dosr():
    preset_path = Path(BASE_DIR, "data", "DosR.txt")
    return preset_path.read_text()


if __name__ == "__main__":
    port = int(os.getenv("PORT", 8080))
    uvicorn.run(app, port=port, log_level="info", host="0.0.0.0", access_log=False)
