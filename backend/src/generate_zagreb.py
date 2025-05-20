import json

from main import generate_map

# Load card configuration
with open("card.json", "r") as f:
    card_config = json.load(f)

# Extract parameters
lat, lon = card_config["center"]
width = card_config["width"]
height = card_config["height"]
location = card_config["location"]
layers = card_config["layers"]

# Set active flag for each layer
for layer in layers:
    layer["active"] = True

# Generate map
formats = ["svg"]  # Output formats
style = "default"  # Map style

generate_map(lat, lon, width, height, formats, style, location, layers)
