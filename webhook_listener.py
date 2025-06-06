from fastapi import FastAPI, Request
# from pydantic import BaseModel
# import uvicorn
import os
import json
import httpx
from protein_ligand_data import make_valid_fname

app = FastAPI()
MAX_RETRIES = 3  # Maximum number of retries for the webhook

# Where to save incoming results
OUTPUT_DIR = "./docking_results"
os.makedirs(OUTPUT_DIR, exist_ok = True)

def write_results(filename, payload):
    filepath = os.path.join(OUTPUT_DIR, filename)
    with open(filepath, "w") as f:
        json.dump(payload, f, indent = 2)

    return filepath

@app.post("/webhook")
async def receive_docking_result(request: Request):
    payload = await request.json()

    status = payload.get("status", "unknown")
    error_type = payload.get("error_type", "")
    result = payload.get("result", {})
    uniprot_id = result.get("uniprot_id", "unknown")
    ligand = result.get("ligand", "unknown")

    retry_count = payload.get("retry_count", 0)

    # If docking failed due to no pose and retry count is less than max retries, resubmit the job
    if status == "failed":
        callback_url = os.getenv("WEBHOOK", "http://your.host") + "/webhook"
        if error_type == "no_pose_generated" and retry_count < MAX_RETRIES:
            print(f"Retrying docking for {uniprot_id} + {ligand} (retry {retry_count + 1})")
            try:
                await httpx.AsyncClient().post(
                    "https://diffdock.toxindex.com/start_docking_uniprot",
                    json = {
                        "uniprot_id": uniprot_id,
                        "ligand": ligand,
                        "callback_url": callback_url,
                        "retry_count": retry_count + 1
                    }
                )
            except Exception as e:
                print(f"Retry submission failed: {e}")

            return {"message": f"Retried docking for {uniprot_id} + {ligand} (retry {retry_count + 1})"}
        elif error_type == "inference_error" and ("CUDA out of memory" in payload.get("stderr", "")):
            if retry_count < MAX_RETRIES:
                try:
                    await httpx.AsyncClient().post(
                        "https://diffdock.toxindex.com/start_docking_uniprot",
                        json = {
                            "uniprot_id": uniprot_id,
                            "ligand": ligand,
                            "callback_url": callback_url,
                            "retry_count": retry_count + 1,
                            "batch_size": 1
                        }
                    )
                except Exception as e:
                    print(f"Retry submission failed: {e}")

                return {"message": f"CUDA out of memory error for {uniprot_id} + {ligand}. Retrying with batch_size = 1."}
                
        else:
            filename = "unknown_unknown.json"
            write_results(filename, payload)
            return {"message": f"Docking failed for {uniprot_id} + {ligand}. No retry."}

    else:
        # Create a unique filename based on UniProt and ligand (you can modify this)
        filename = make_valid_fname(uniprot_id, ligand)
        filepath = write_results(filename, payload)

        return {"message": f"Result saved to {filepath}"}
