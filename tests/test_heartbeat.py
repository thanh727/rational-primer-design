import time
from fastapi.testclient import TestClient
import rational_design.api as api_module

def test_heartbeat_endpoint():
    app = api_module.create_app()
    client = TestClient(app)

    # Initial state
    assert api_module.HAS_HEARTBEAT_STARTED is False

    # Call heartbeat endpoint
    response = client.post("/api/heartbeat")
    assert response.status_code == 200
    assert response.json() == {"status": "alive"}

    # Assert global state updated
    assert api_module.HAS_HEARTBEAT_STARTED is True
    assert api_module.LAST_HEARTBEAT <= time.time()
