import os
import uvicorn
from fastapi import FastAPI, status
import config.base_models as base_models
from manager import Manager


app = FastAPI()


@app.get("/", status_code=status.HTTP_200_OK)
def root():
    return {"INFO": "This is the root of motif-search service for finding motif in given DNA sequences"}


@app.post("/find-motif")
async def find_motif(request_data: base_models.FindMofif):
    manager = Manager(request_data)
    response = manager.find_motif()
    return response


if __name__ == "__main__":
    port = int(os.getenv("PORT", 8080))
    uvicorn.run(app, port=port, log_level="info", host="0.0.0.0", access_log=False)
