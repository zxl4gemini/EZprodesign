"""Launch the local-only EZprodesign preview server."""
from . import create_app
if __name__ == "__main__":
    create_app().run(host="127.0.0.1", port=18785, debug=False, use_reloader=False)
