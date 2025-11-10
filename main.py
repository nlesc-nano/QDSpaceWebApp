import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from mangum import Mangum
from backend.app import app

handler = Mangum(app, api_gateway_base_path="/prod/")
logger.info("Handler initialized with dspath=/prod/")