"""Launch the loopback-only EZprodesign preview server."""
import argparse
from . import create_app

def main():
    parser = argparse.ArgumentParser(description="Run the local EZprodesign browser workbench.")
    parser.add_argument("--port", type=int, default=18785, help="loopback port (default: 18785)")
    args = parser.parse_args()
    if not 1024 <= args.port <= 65535:
        parser.error("port must be between 1024 and 65535")
    create_app().run(host="127.0.0.1", port=args.port, debug=False, use_reloader=False)

if __name__ == "__main__":
    main()
