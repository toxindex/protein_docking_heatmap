from fastapi import FastAPI, Request
# from pydantic import BaseModel
# import uvicorn
import os
import json
from protein_ligand_data import make_valid_fname

app = FastAPI()

# Where to save incoming results
OUTPUT_DIR = "./docking_results"
os.makedirs(OUTPUT_DIR, exist_ok = True)

@app.post("/webhook")
async def receive_docking_result(request: Request):
    payload = await request.json()

    status = payload.get("status", "unknown")
    result = payload.get("result", {})
    uniprot_id = result.get("uniprot_id", "unknown")
    ligand = result.get("ligand", "unknown")

    # Create a unique filename based on UniProt and ligand (you can modify this)
    filename = make_valid_fname(uniprot_id, ligand)
    filepath = os.path.join(OUTPUT_DIR, filename)

    with open(filepath, "w") as f:
        json.dump(payload, f, indent = 2)

    return {"message": f"Result saved to {filepath}"}
