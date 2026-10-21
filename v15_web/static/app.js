'use strict';
(() => {
const $ = id => document.getElementById(id);
const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt = (x,n=2) => Number.isFinite(x) ? x.toFixed(n) : '未知';
const flags = {
 EXPLORATORY_GEOMETRY_NOT_BINDING_OR_AFFINITY:'仅探索性几何分析，未验证结合或亲和力',
 GLYCAN_OCCUPANCY_UNKNOWN:'糖基化占有率未知，序列基序不等于实际糖链',
 MEMBRANE_ORIENTATION_UNKNOWN:'膜朝向未知', MEMBRANE_ACCESS_NOT_CERTIFIED:'膜环境可及性未认证',
 UNVERIFIED_CONFIDENCE:'结构置信度未核验', PLDDT_UNVERIFIED:'无已核验的 pLDDT，未将 B 因子冒充置信度',
 COARSE_ENVELOPE_NOT_DOCKED_ANTIBODY:'抗体包络仅为几何近似，尚未真实对接',
 MEMBRANE_MICELLE_ACCESS_UNCERTIFIED:'膜 / 胶束环境下的接近通道未认证',
 NO_COHERENT_COLLISION_FREE_APPROACH:'未找到共同无碰撞接近方向',
 HIGH_HYDROPHOBIC_CONTEXT_RISK:'高疏水表面，可能被胶束或膜遮挡',
 HYDROPHOBIC_CONTEXT_RISK_UNCERTIFIED:'局部疏水风险，真实遮挡情况未知',
 POTENTIAL_GLYCAN_NEAR_HOTSPOTS_OCCUPANCY_UNKNOWN:'热点邻近潜在糖基化位点，实际占有率未知',
 GLYCO_SEQUON_ASN_OCCUPANCY_UNKNOWN_EXCLUDED_HOTSPOT:'糖基化基序 Asn 已排除为热点，糖链占有率未知',
 NEAR_POTENTIAL_GLYCAN:'邻近潜在糖链', INCOMPLETE_BACKBONE_OR_SIDECHAIN:'主链或侧链原子不完整',
 FEWER_THAN_THREE_SUPPORTED_PATCHES_NO_FORCED_RECOMMENDATIONS:'证据支持的分组少于三个，不补造热点'
};
const flagHTML = values => [...new Set(values || [])].map(v => `<div class="flag">${esc(flags[v] || v)}<br><code>${esc(v)}</code></div>`).join('');
let file=null, result=null, viewer=null, model=null, surface=null, surfaceGeneration=0, revision=0, controller=null, selected=null, applied=new Set(), spinning=false, busy=false, renderError=null, renderReady=false, maxFile=15*1024*1024;
const records = () => new Map((result?.residues || []).map(r => [r.resi,r]));
const mapping = () => new Map((result?.target?.numbering || []).map(r => [r.resi,r]));
const originalLabel = n => {const m=mapping().get(n); return m ? `${m.chain_original || '∅'}:${m.resseq}${m.icode || ''}` : `A:${n}`;};
function color(atom){const r=recordsCache.get(atom.resi);const x=Math.max(0,Math.min(1,(r?.color_score || 0)/.85));const stops=[[200,81,73],[233,214,189],[169,203,212],[40,95,158]];const u=x*3,i=Math.min(2,Math.floor(u)),f=u-i;const c=stops[i].map((v,j)=>Math.round(v+(stops[i+1][j]-v)*f));return (c[0]<<16)+(c[1]<<8)+c[2];}
let recordsCache=new Map();
function clearScene(){surfaceGeneration++;renderReady=false;if(viewer){viewer.spin(false);viewer.removeAllSurfaces();viewer.removeAllModels();viewer.removeAllShapes();viewer.removeAllLabels();viewer.render();}surface=null;model=null;spinning=false;$('rotate').textContent='旋转';}
function invalidate(){busy=false;$('progress').classList.remove('busy');result=null;selected=null;applied.clear();clearScene();$('result').hidden=true;document.body.classList.remove('has-result');updateExports();}
function fail(message){$('error').textContent=message;$('error').hidden=false;}
function populateChains(chains, value='auto'){$('targetChain').replaceChildren(new Option('自动选择','auto'),...chains.map(c=>new Option(`链 ${c || '∅'}`,c)));$('targetChain').value=value;}
async function chooseFile(next){
 if(!next)return;
 controller?.abort();revision++;invalidate();file=null;$('analyze').disabled=true;$('error').hidden=true;populateChains([]);
 if(!/\.(pdb|cif|mmcif|fa|fasta|faa)$/i.test(next.name)){fail('支持 PDB、mmCIF / AF3 mmCIF 或 FASTA；请解压后选择结构文件。');return;}
 if(!next.size || next.size>maxFile){fail(`文件为空或超过引擎 ${Math.floor(maxFile/1024/1024)} MiB 上限。`);return;}
 file=next;$('filename').textContent=next.name;$('analyze').disabled=false;await analyze();
}
async function analyze(){
 if(!file)return;
 controller?.abort();controller=new AbortController();const token=++revision;invalidate();busy=true;$('error').hidden=true;$('progress').classList.add('busy');$('progress').textContent='正在本机解析结构、计算可及表面和共同接近通道…';$('analyze').disabled=true;
 const body=new FormData();body.append('file',file);body.append('environment',$('environment').value);body.append('binder_type',$('binderType').value);body.append('target_chain',$('targetChain').value);
 try{const response=await fetch('/api/analyze',{method:'POST',body,signal:controller.signal});const data=await response.json();if(token!==revision)return;if(!response.ok)throw new Error(data.error || `分析失败 (${response.status})`);result=data;recordsCache=records();applied=new Set(data.patches.filter(p=>p.status==='DESIGN_EXPLORATION').map(p=>p.id));populateChains(data.source.parsed_chain_names || [],data.provenance.parameters.target_chain);await showResult();if(token!==revision)return;$('progress').textContent=data.status==='needs_structure'?'序列已读取；需要三维结构后才能定位热点。':'分析完成。可点击热点查看实际残基，勾选分组后导出复核清单。';}
 catch(error){if(token!==revision || error.name==='AbortError')return;invalidate();fail(error.message);$('progress').textContent='分析未完成，请检查文件或参数后重试。';}
 finally{if(token===revision){busy=false;$('progress').classList.remove('busy');$('analyze').disabled=!file;}}
}
async function showResult(){
 const d=result,t=d.target;$('result').hidden=false;document.body.classList.add('has-result');$('resultTitle').textContent=d.source.filename;$('resultMeta').textContent=`${t.residue_count} aa · ${t.atom_count} 原子 · 原始链 ${t.chain_original ?? '无结构链'} · ${d.source.format.toUpperCase()} · ${d.provenance.parameters.binder_type}`;
 const explore=d.patches.filter(p=>p.status==='DESIGN_EXPLORATION').length;$('statusBadge').textContent=d.status==='needs_structure'?'三维结构待补充':`${explore} 个探索分组 · 未验证结合`;
 $('environmentNote').textContent={detergent:'DDM / GDN 胶束：动态包覆层的位置与形状未知。仅根据暴露疏水性提示遮挡风险，不绘制精确胶束边界或膜平面。',soluble:'可溶环境：在上传结构背景中评估目标链表面。未上传的复合物组分、真实糖链和构象变化仍可能改变可及性。',bilayer:'细胞双层膜：朝向与插膜范围尚未注册。不推断固定膜平面；所有分组必须保留待复核状态。'}[d.environment.mode];
 $('patchCount').textContent=`${d.patches.length} 组 / ${explore} 探索`;$('patches').innerHTML=d.patches.map(p=>`<article class="patch-card" data-patch="${esc(p.id)}"><button class="focus-patch" data-focus="${esc(p.id)}" aria-label="聚焦 ${esc(p.id)}"></button><div class="patch-top"><input type="checkbox" data-apply="${esc(p.id)}" aria-label="导出 ${esc(p.id)}" ${applied.has(p.id)?'checked':''}><strong>${esc(p.id)} · ${p.members.length} 个残基</strong><span class="score">${fmt(p.score*100,1)}</span></div><div class="patch-status ${p.status==='HOLD_REVIEW'?'hold':''}">${p.status==='HOLD_REVIEW'?'待复核 HOLD_REVIEW':'设计探索 DESIGN_EXPLORATION'}</div><div class="residue-chips">${p.hotspots.map(n=>`${esc(recordsCache.get(n)?.aa)} ${esc(originalLabel(n))}`).join(' · ')}</div><small>热点 / 原始编号 · 分数为启发式优先级</small></article>`).join('');
 $('emptyPatches').hidden=!!d.patches.length;$('emptyPatches').textContent=d.status==='needs_structure'?'FASTA 仅有序列。未运行结构预测；没有三维坐标，因此不生成热点。请上传实验结构或预测 PDB / mmCIF。':`没有支持的成组热点。暴露残基 ${d.qc.exposed_residue_count ?? '未知'} 个，通过单残基接近探针 ${d.qc.approach_clear_residue_count ?? '未知'} 个。原子不完整 ${d.residues.filter(r=>(r.warnings||[]).includes('INCOMPLETE_BACKBONE_OR_SIDECHAIN')).length} 个；低暴露 ${d.residues.filter(r=>r.access_status==='LOW_EXPOSURE').length} 个；接近受阻 ${d.residues.filter(r=>r.access_status==='APPROACH_BLOCKED').length} 个。分组还须满足至少两个可用热点与连续表面；不强制补造。`;
 $('patchDetail').textContent=d.patches.length?'选择一个热点分组查看证据。':$('emptyPatches').textContent;$('residueDetail').textContent='';$('selectedLabel').textContent='';
 $('warnings').innerHTML=flagHTML(d.warnings);$('qc').textContent=JSON.stringify(d.qc,null,2);$('sequence').textContent=t.sequence;$('sequenceNote').textContent=(t.glyco_sequons || []).length ? `潜在 N-糖基化序列基序（不代表占有率）：${t.glyco_sequons.map(x=>`${x.resi} ${x.motif}`).join('、')}`:'序列基序与真实糖基化、折叠及结合能力不能等同。';$('sourceHash').textContent=`源文件 SHA-256: ${d.source.sha256} | 引擎 ${d.provenance.engine_version}: ${d.provenance.engine_sha256}`;
 if(t.pdb){try{await draw();}catch(e){renderError=String(e);$('viewerEmpty').hidden=false;$('viewerEmpty').textContent='3D 渲染不可用；分析证据与导出仍可查看。';fail(`3D 渲染失败：${e.message}`);}}else{$('viewerEmpty').hidden=false;$('viewerEmpty').textContent='FASTA 序列已读取 · 三维结构待补充';}
 updateExports();
}
async function draw(){
 const token=revision;renderError=null;if(!window.$3Dmol)throw new Error('本地 3Dmol 未加载');if(!viewer)viewer=$3Dmol.createViewer($('viewer'),{backgroundColor:'#f4f7f5',antialias:true});
 model=viewer.addModel(result.target.pdb,'pdb');model.setStyle({},{cartoon:{colorfunc:color}});model.setClickable({},true,atom=>inspectResidue(atom.resi));viewer.zoomTo();viewer.render();viewer.resize();$('viewerEmpty').hidden=true;await drawSurface();if(token!==revision)return;renderReady=true;if(result?.patches.length)focusPatch(result.patches[0].id,false);
}
async function drawSurface(){
 const generation=++surfaceGeneration;if(surface!==null&&viewer){viewer.removeSurface(surface);surface=null;}if(!viewer||!model||!$('surfaceToggle').checked){viewer?.render();return;}
 // Synchronous CPU surface avoids worker/CDN fetches; capped by engine atom limit.
 $3Dmol.setSyncSurface(true);
 const added=viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:.63,colorfunc:color},{},{},{});surface=added.surfid;await added;
 if(generation!==surfaceGeneration)return;viewer.render();
}
function focusPatch(id,zoom=true){
 const p=result?.patches.find(x=>x.id===id);if(!p||!viewer||!model)return;selected=id;
 document.querySelectorAll('.patch-card').forEach(e=>e.classList.toggle('active',e.dataset.patch===id));
 model.setStyle({},{cartoon:{colorfunc:color}});model.setStyle({chain:'A',resi:p.members},{stick:{radius:.14,color:'#48867a'}},true);model.setStyle({chain:'A',resi:p.hotspots},{stick:{radius:.25,color:'#e8b24c'},sphere:{scale:.26,color:'#e8b24c'}},true);
 viewer.removeAllLabels();viewer.removeAllShapes();for(const n of p.hotspots)viewer.addLabel(`${recordsCache.get(n)?.aa}${originalLabel(n)}`,{fontSize:11,backgroundColor:'#15352f',backgroundOpacity:.8,fontColor:'white',inFront:true},{chain:'A',resi:n,atom:'CA'});
 if(p.approach_direction&&p.center){const start={x:p.center[0],y:p.center[1],z:p.center[2]},end={x:start.x+p.approach_direction[0]*15,y:start.y+p.approach_direction[1]*15,z:start.z+p.approach_direction[2]*15};viewer.addArrow({start,end,radius:.25,color:'#368579'});}
 if(zoom){viewer.spin(false);spinning=false;$('rotate').textContent='旋转';viewer.zoomTo({chain:'A',resi:p.members});
 const direction=p.approach_direction;if(direction){const len=Math.hypot(...direction);if(len>0){const [x,y,z]=direction.map(v=>v/len);let q=z<-.999999?[1,0,0,0]:[y,-x,0,1+z];const norm=Math.hypot(...q);q=q.map(v=>v/norm);const view=viewer.getView();q.forEach((v,i)=>view[i+4]=v);viewer.setView(view);}}} viewer.render();$('selectedLabel').textContent=id;
 const rr=p.members.map(n=>recordsCache.get(n)).filter(Boolean),area=rr.reduce((a,r)=>a+(r.sidechain_sas_A2||0),0),risk=rr.reduce((a,r)=>a+(r.hydrophobic_risk||0),0)/Math.max(1,rr.length),glyco=rr.filter(r=>r.glyco_sequon).length;
 const physical=p.access_checks?.physical_pass?'通过粗粒度检查':'未通过';
 $('patchDetail').innerHTML=`<div class="metrics"><div class="metric">组内侧链 SASA 总和<b>${fmt(area,1)} Å²</b></div><div class="metric">平均疏水风险代理<b>${fmt(risk,2)}</b></div><div class="metric">糖基化基序 Asn<b>${glyco} 个</b></div><div class="metric">共同接近通道<b>${physical}</b></div><div class="metric">最小采样间隙<b>${fmt(p.access_checks?.sampled_min_clearance_A)} Å</b></div><div class="metric">结构置信度<b>未核验</b></div></div><p>箭头仅表示引擎检查的接近方向。侧链面积为组内求和，不是抗体埋藏界面面积。</p><p>中心 / Å：${(p.center||[]).map(x=>fmt(x,3)).join(', ')} · ${esc(p.binder_type)} · ${esc(p.status)}</p><div class="table-scroll"><table><thead><tr><th>原始链:编号</th><th>归一化</th><th>残基</th><th>侧链 Å²</th><th>RSA 代理</th><th>Cα 坐标 / Å</th></tr></thead><tbody>${rr.map(r=>`<tr><td><button data-residue="${r.resi}">${esc(originalLabel(r.resi))}${p.hotspots.includes(r.resi)?' ★':''}</button></td><td>A:${r.resi}</td><td>${esc(r.aa)}</td><td>${fmt(r.sidechain_sas_A2,1)}</td><td>${fmt(r.rsa_proxy)}</td><td>${(r.ca||[]).map(x=>fmt(x,2)).join(', ')}</td></tr>`).join('')}</tbody></table></div>${flagHTML(p.warnings)}<details><summary>原始几何检查与分数组成</summary><pre>${esc(JSON.stringify({rationale:p.rationale,access_checks:p.access_checks,approach_direction:p.approach_direction},null,2))}</pre></details>`;
}
function inspectResidue(n){const r=recordsCache.get(n);if(!r)return;$('residueDetail').innerHTML=`<p><b>${esc(r.aa)} ${esc(originalLabel(n))} → A:${n}</b> · ${esc(r.access_status)} · ${r.accessible_probe_count} 个通过的接近探针</p><p>Cα / Å：${(r.ca||[]).map(x=>fmt(x,3)).join(', ')} · pLDDT：${r.plddt_verified ? fmt(r.plddt,1) : '未核验'}</p>${flagHTML(r.warnings)}`;viewer?.zoomTo({chain:'A',resi:[n]},300);viewer?.render();}
function updateExports(){const valid=!!result;const count=applied.size;$('exportInfo').textContent=valid?`已勾选 ${count} / ${result.patches.length} 个分组。所有清单均为人工复核用途。`:'请选择需要保留的分组。';$('exportJson').disabled=!valid;$('exportTsv').disabled=!valid;$('exportPml').disabled=!result?.target?.pdb;}
function manifest(){if(!result)throw new Error('尚无当前结果');return {schema:'EZprodesign.v1.5.review',status:result.status,handoff_type:'HUMAN_REVIEW_DESIGN_EXPLORATION',binding_validated:false,affinity_validated:false,applied_patch_ids:result.patches.filter(p=>applied.has(p.id)).map(p=>p.id),source:result.source,provenance:result.provenance,environment:result.environment,target:result.target,patches:result.patches.filter(p=>applied.has(p.id)).map(p=>({...p,member_mapping:p.members.map(n=>({...mapping().get(n),...recordsCache.get(n),is_hotspot:p.hotspots.includes(n)}))})),warnings:result.warnings,qc:result.qc,uncertainty:'Unverified structure confidence and glycan occupancy. Geometry proxies are not docking, affinity, micelle boundary, or certified membrane access.'};}
function save(text,name,type){const url=URL.createObjectURL(new Blob([text],{type}));const a=document.createElement('a');a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),3000);}
const safeName = () => (result.source.filename.replace(/\.[^.]*$/,'').replace(/[^a-zA-Z0-9_-]/g,'_').slice(0,70)||'target');
function exportTSV(){const m=manifest();const cols=['patch_id','status','is_hotspot','normalized_chain','normalized_resi','original_chain','original_resseq','original_icode','aa','ca_x_A','ca_y_A','ca_z_A','sidechain_sas_A2','rsa_proxy','hydrophobic_risk','plddt_verified','approach_direction','access_checks','warnings','source_sha256','engine_version','engine_sha256','environment','binder_type'];const clean=x=>{const s=typeof x==='object'?JSON.stringify(x):String(x??'');return (/^[=+@-]/.test(s)&&!/^[-+]?\d/.test(s)?"'"+s:s).replace(/[\t\r\n]/g,' ');};const rows=m.patches.flatMap(p=>p.member_mapping.map(r=>[p.id,p.status,r.is_hotspot,'A',r.resi,r.chain_original,r.resseq,r.icode,r.aa,...(r.ca||['','','']),r.sidechain_sas_A2,r.rsa_proxy,r.hydrophobic_risk,r.plddt_verified,p.approach_direction,p.access_checks,[...m.warnings,...p.warnings,...r.warnings],m.source.sha256,m.provenance.engine_version,m.provenance.engine_sha256,m.environment.mode,p.binder_type]));save([cols,...rows].map(row=>row.map(clean).join('\t')).join('\n')+'\n',`EZv15_${safeName()}_review.tsv`,'text/tab-separated-values');}
function base64(s){const bytes=new TextEncoder().encode(s);let out='';for(let i=0;i<bytes.length;i+=8192)out+=String.fromCharCode(...bytes.subarray(i,i+8192));return btoa(out);}
function exportPML(){const m=manifest();const payload=base64(JSON.stringify(m));const code=`# EZprodesign v1.5 human review only\n# Source ${safeName()}\n# Geometry only. Binding and affinity are not validated.\npython\nimport base64, json, re\nfrom pymol import cmd\nd = json.loads(base64.b64decode('${payload}').decode('utf-8'))\ncmd.read_pdbstr(d['target']['pdb'], 'ezv15_target')\ncmd.hide('everything', 'ezv15_target')\ncmd.show('cartoon', 'ezv15_target')\ncmd.color('gray70', 'ezv15_target')\ncmd.show('surface', 'ezv15_target')\ncmd.set('transparency', 0.55, 'ezv15_target')\nfor index, p in enumerate(d['patches'], 1):\n    members = '+'.join(str(int(n)) for n in p['members'])\n    hot = '+'.join(str(int(n)) for n in p['hotspots'])\n    safe_id = re.sub('[^a-zA-Z0-9_]', '_', p['id'])\n    name = 'ezv15_patch_' + safe_id\n    count = cmd.select(name, 'ezv15_target and chain A and resi ' + members)\n    cmd.show('sticks', name)\n    cmd.color('marine', name)\n    hname = 'ezv15_hotspots_' + safe_id\n    hcount = cmd.select(hname, 'ezv15_target and chain A and resi ' + hot)\n    cmd.color('yelloworange', hname)\n    print('EZV15_PATCH', p['id'], count, hcount)\ncmd.orient('ezv15_target')\nprint('EZV15_TARGET_ATOMS', cmd.count_atoms('ezv15_target'))\nprint('EZV15_APPLIED_IDS', json.dumps(d['applied_patch_ids']))\npython end\n`;save(code,`EZv15_${safeName()}_review.pml`,'text/plain');}
$('fileInput').addEventListener('change',e=>chooseFile(e.target.files[0]));for(const name of ['dragenter','dragover'])$('dropzone').addEventListener(name,e=>{e.preventDefault();$('dropzone').classList.add('dragover');});for(const name of ['dragleave','drop'])$('dropzone').addEventListener(name,e=>{e.preventDefault();$('dropzone').classList.remove('dragover');if(name==='drop')chooseFile(e.dataTransfer.files[0]);});
$('analyze').addEventListener('click',analyze);for(const id of ['environment','binderType','targetChain'])$(id).addEventListener('change',()=>{if(file)analyze();});
async function sample(name){const token=++revision;controller?.abort();invalidate();file=null;$('analyze').disabled=true;$('error').hidden=true;try{const r=await fetch(`/api/sample/${name}`);if(!r.ok)throw new Error('此本地示例尚未就绪。');const blob=await r.blob();if(token===revision)await chooseFile(new File([blob],`${name}.pdb`,{type:'chemical/x-pdb'}));}catch(e){if(token===revision)fail(e.message);}}
$('sampleUbq').addEventListener('click',()=>sample('1ubq'));
$('patches').addEventListener('click',e=>{const el=e.target.closest('[data-focus]');if(el)focusPatch(el.dataset.focus);});$('patches').addEventListener('change',e=>{const id=e.target.dataset.apply;if(!id)return;e.target.checked?applied.add(id):applied.delete(id);updateExports();});$('patchDetail').addEventListener('click',e=>{const el=e.target.closest('[data-residue]');if(el)inspectResidue(Number(el.dataset.residue));});
$('surfaceToggle').addEventListener('change',()=>drawSurface().catch(e=>fail(e.message)));$('rotate').addEventListener('click',()=>{spinning=!spinning;viewer?.spin(spinning?'y':false,.5);$('rotate').textContent=spinning?'停止旋转':'旋转';});$('resetView').addEventListener('click',()=>{viewer?.zoomTo({},400);viewer?.render();});window.addEventListener('resize',()=>viewer?.resize());
$('exportJson').addEventListener('click',()=>save(JSON.stringify(manifest(),null,2)+'\n',`EZv15_${safeName()}_review.json`,'application/json'));$('exportTsv').addEventListener('click',exportTSV);$('exportPml').addEventListener('click',exportPML);
fetch('/api/health').then(r=>r.json()).then(d=>{$('serviceState').textContent=d.engine_available?'本地分析引擎已就绪':'等待本地分析引擎';$('sampleUbq').disabled=!d.samples['1ubq'];if(d.max_file_bytes){maxFile=d.max_file_bytes;$('fileLimit').textContent=`本机内存处理 · 文件上限 ${Math.floor(maxFile/1024/1024)} MiB / 请求上限 16 MiB · 不自动预测结构`;}}).catch(()=>{$('serviceState').textContent='本地服务不可达';});
function webglState(){const gl=viewer?.getRenderer()?.getContext();return gl?{version:gl.getParameter(gl.VERSION),lost:gl.isContextLost(),width:gl.drawingBufferWidth,height:gl.drawingBufferHeight}:null;}
// Read-only inspection hooks for reproducible browser checks; no persistent storage.
window.EZV15={getResult:()=>result,getManifest:()=>result?manifest():null,getState:()=>({busy,revision,selected,applied:[...applied],renderReady,renderError,webgl:webglState(),atoms:model?.selectedAtoms({}).length||0,highlightedAtoms:model?.selectedAtoms({}).filter(a=>a.style?.sphere?.color==='#e8b24c').map(a=>({resi:a.resi,x:a.x,y:a.y,z:a.z}))||[],view:viewer?.getView()}),getAtoms:()=>model?.selectedAtoms({}).map(a=>({resi:a.resi,atom:a.atom,chain:a.chain,x:a.x,y:a.y,z:a.z}))||[]};
})();
