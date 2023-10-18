import os

current_route = os.path.abspath(__file__)
parent_route = os.path.dirname(current_route)
root_route = os.path.dirname(parent_route)
save_route_images = os.path.join(root_route, "images")
save_route_data = os.path.join(root_route, "data")
