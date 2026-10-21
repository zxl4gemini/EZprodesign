"""Local-only, in-memory upload web application for EZprodesign 1.5."""
from __future__ import annotations
import hashlib
import io
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlsplit
from flask import Flask, Request, jsonify, render_template, request, send_file
from werkzeug.exceptions import RequestEntityTooLarge

ROOT = Path(__file__).resolve().parents[1]
MAX_REQUEST = 16 * 1024 * 1024

class MemoryRequest(Request):
    # Werkzeug's default spools large multipart files to disk. Never do that here.
    def _get_file_stream(self, total_content_length, content_type, filename=None, content_length=None):
        return io.BytesIO()


@lru_cache(maxsize=1)
def load_engine():
    from v15_engine import core
    return core, hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest()

def sample_paths():
    # Only a publicly licensed example; never load private laboratory data.
    fixture = Path(__file__).resolve().parent / 'examples/1UBQ.pdb'
    return {'1ubq': fixture if fixture.is_file() else None}

def create_app():
    app = Flask(__name__)
    app.request_class = MemoryRequest
    app.config.update(MAX_CONTENT_LENGTH=MAX_REQUEST, MAX_FORM_MEMORY_SIZE=MAX_REQUEST,
                      MAX_FORM_PARTS=12, JSON_SORT_KEYS=False)
    app.json.ensure_ascii = False

    @app.before_request
    def local_only():
        if request.host.split(':')[0] not in ('localhost', '127.0.0.1'):
            return jsonify(error='只允许 localhost / 127.0.0.1 本机访问。'), 403
        if request.method == 'POST':
            origin = request.headers.get('Origin')
            if origin and urlsplit(origin).netloc != request.host:
                return jsonify(error='拒绝跨站上传，请从本机页面重试。'), 403

    @app.after_request
    def headers(response):
        response.headers['Cache-Control'] = 'no-store'
        response.headers['X-Content-Type-Options'] = 'nosniff'
        response.headers['Referrer-Policy'] = 'no-referrer'
        response.headers['Content-Security-Policy'] = "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; worker-src 'self' blob:; object-src 'none'; base-uri 'none'; frame-ancestors 'none'"
        return response

    @app.errorhandler(RequestEntityTooLarge)
    def too_large(_):
        return jsonify(error='上传请求超过 16 MiB。请减小文件后重试。'), 413

    @app.get('/')
    def index():
        return render_template('index.html')

    @app.get('/vendor/3Dmol-min.js')
    def vendor():
        return send_file(Path(__file__).resolve().parent / 'vendor/3Dmol-min.js', mimetype='text/javascript')

    @app.get('/api/health')
    def health():
        try:
            core, engine_hash = load_engine()
            available, limit = True, core.MAX_BYTES
        except ImportError:
            available, limit, engine_hash = False, None, None
        return jsonify(status='ok' if available else 'engine_pending', app_version='1.5.0',
                       engine_available=available, engine_sha256=engine_hash, max_request_bytes=MAX_REQUEST, max_file_bytes=limit,
                       upload_persistence=False, samples={k: bool(v) for k,v in sample_paths().items()})

    @app.get('/api/sample/<sample>')
    def sample(sample):
        path = sample_paths().get(sample)
        if path is None:
            return jsonify(error='本机没有此示例结构，您可以上传自己的文件。'), 404
        return send_file(path, mimetype='chemical/x-pdb', download_name=sample+'.pdb')

    @app.post('/api/analyze')
    def analyze():
        if 'file' not in request.files:
            return jsonify(error='请选择 PDB、mmCIF、AF3 mmCIF 或 FASTA 文件。'), 400
        upload = request.files['file']
        filename = (upload.filename or '').replace('\\', '/').rsplit('/', 1)[-1]
        if not filename or len(filename) > 240:
            return jsonify(error='文件名为空或过长。'), 400
        environment = request.form.get('environment', 'detergent')
        binder_type = request.form.get('binder_type', 'VHH')
        target_chain = request.form.get('target_chain', 'auto')
        if environment not in ('detergent','soluble','bilayer') or binder_type not in ('VHH','Fab'):
            return jsonify(error='环境或抗体类型无效。'), 400
        if not target_chain or len(target_chain) > 32 or any(ord(c) < 32 for c in target_chain):
            return jsonify(error='目标链标识无效。'), 400
        try:
            core, engine_hash = load_engine()
        except ImportError:
            return jsonify(error='分析引擎尚未就绪，请稍后重试。'), 503
        raw = upload.read(MAX_REQUEST + 1)
        if len(raw) > core.MAX_BYTES:
            return jsonify(error=f'分析引擎文件上限为 {core.MAX_BYTES // (1024*1024)} MiB，请减小文件。'), 413
        try:
            result = core.analyze_bytes(filename, raw, target_chain=target_chain,
                                        environment=environment, binder_type=binder_type)
        except ValueError as exc:
            return jsonify(error='无法分析此文件：' + str(exc)), 400
        except Exception:
            # Never log upload bytes or detailed parsing exceptions to disk.
            return jsonify(error='本机分析未完成。请检查结构完整性或尝试较小的单链结构。'), 500
        result['provenance'] = {
            'engine_version': result.get('qc',{}).get('engine_version', 'unknown'),
            'engine_sha256': engine_hash,
            'web_version': '1.5.0', 'upload_persisted': False,
            'parameters': {'environment': environment, 'binder_type': binder_type, 'target_chain': target_chain},
        }
        return jsonify(result)

    return app
