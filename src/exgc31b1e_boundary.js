function run31E(){
 const seeds={single:[[0,1]],path:[[0,1],[1,2],[2,3]],triangle:[[0,1],[1,2],[2,0]],cycle:[[0,1],[1,2],[2,3],[3,0]],star:[[0,1],[0,2],[0,3]],K4:[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]};
 const pair=(a,b)=>a<b?"("+a+"|"+b+")":"("+b+"|"+a+")";
 const cmp=(a,b)=>a<b?-1:a>b?1:0;
 function mincut(V,B,roots){
  let A=[...V.values()],ids=new Map(A.map((v,i)=>[v.id,i])),N=A.length,src=N,sink=N+1,adj=Array.from({length:N+2},()=>[]),outside=new Set();
  function link(a,b,c){let x={to:b,cap:c,rev:adj[b].length},y={to:a,cap:0,rev:adj[a].length};adj[a].push(x);adj[b].push(y)}
  for(let b of B.values()){let [a,c]=b.ends;let i=ids.get(a),j=ids.get(c);link(i,j,1);link(j,i,1);if(A[i].root!==A[j].root){outside.add(a);outside.add(c)}}
  let left=A.filter(v=>v.root===roots[0]&&!outside.has(v.id)),right=A.filter(v=>v.root===roots[1]&&!outside.has(v.id));
  if(!left.length||!right.length)return {flow:null,left:left.length,right:right.length};
  for(let v of left)link(src,ids.get(v.id),B.size+1);for(let v of right)link(ids.get(v.id),sink,B.size+1);
  let flow=0;while(true){let prev=Array(N+2).fill(null),q=[src];prev[src]=[-1,-1];for(let j=0;j<q.length&&prev[sink]===null;j++){let a=q[j];for(let k=0;k<adj[a].length;k++){let e=adj[a][k];if(e.cap>0&&prev[e.to]===null){prev[e.to]=[a,k];q.push(e.to);}}}if(prev[sink]===null)break;let d=Infinity;for(let b=sink;b!==src;){let[a,k]=prev[b];d=Math.min(d,adj[a][k].cap);b=a}for(let b=sink;b!==src;){let[a,k]=prev[b],e=adj[a][k];e.cap-=d;adj[b][e.rev].cap+=d;b=a}flow+=d}
  return {flow,left:left.length,right:right.length};
 }
 function run(seed,R,mode,rs,perm,flowOn){
  let V=new Map,B=new Map,H=new Map,pending=new Map,retiredV=0,retiredB=0,crossJoins=0,maxD=0,maxEnabled=0,adjacent=false,rng=rs>>>0;
  function vKey(v){return "v"+v.id}
  function bKey(b){return "b"+b.id}
  function offerV(v){if(v.g<R&&v.bonds.length===v.need)pending.set(vKey(v),{type:"v",x:v});}
  function vertex(id,n,root,g,birth,need){if(V.has(id))throw Error("duplicate v");if(need<1||need>3)throw Error("degree reservation");maxD=Math.max(maxD,need);let v={id,n,root,g,birth,need,bonds:[]};V.set(id,v);return v;}
  function bond(a,b,g,birth){let id=pair(a.id,b.id),n=pair(a.n,b.n);if(B.has(id))throw Error("duplicate bond");if(a.g!==g||b.g!==g)throw Error("mixed generation");let e={id,n,g,birth,ends:[a.id,b.id],roots:[a.root,b.root],offers:new Map};B.set(id,e);for(let v of [a,b]){v.bonds.push(id);if(v.bonds.length>v.need)throw Error("active degree");offerV(v)}return e}
  let roots=[...new Set(seed.flat())].sort((a,b)=>a-b);
  for(let r of roots)vertex("s"+perm[r],"s"+r,r,0,null,seed.filter(e=>e.includes(r)).length);
  for(let [a,b] of seed)bond(V.get("s"+perm[a]),V.get("s"+perm[b]),0,null);
  while(pending.size){
   maxEnabled=Math.max(maxEnabled,pending.size);
   if(!adjacent)for(let p of pending.values()){if(p.type==="v"&&p.x.bonds.some(id=>B.get(id).ends.some(v=>v!==p.x.id&&pending.has("v"+v)))){adjacent=true;break}}
   let ks=[...pending.keys()].sort(cmp),k=mode==="asc"?ks[0]:mode==="desc"?ks.at(-1):ks[(rng=(Math.imul(1664525,rng)+1013904223)>>>0)%ks.length];let {type,x}=pending.get(k);pending.delete(k);
   if(type==="v"){
    let v=x;if(!V.has(v.id)||v.bonds.length!==v.need)throw Error("bad v");let d=v.need,ng=v.g+1,ev="v"+v.n,bs=v.bonds.map(b=>B.get(b)).sort((a,b)=>cmp(a.id,b.id));
    H.set(k,{id:ev,kind:type,g:v.g,root:v.root,parents:[...new Set([v.birth,...bs.map(b=>b.birth)].filter(Boolean))].sort(),inputs:bs.map(b=>b.n).sort()});retiredV++;V.delete(v.id);
    let D=bs.map(b=>[0,1].map(l=>vertex(v.id+"<"+b.id+":"+l+">",v.n+"<"+b.n+":"+l+">",v.root,ng,ev,d===1?2:3)));
    for(let i=0;i<d;i++)for(let j=0;j<d;j++)if(d<3||i!==j)bond(D[i][0],D[j][1],ng,ev);
    for(let i=0;i<d;i++){let b=bs[i];if(b.offers.has(v.id))throw Error("duplicate offer");b.offers.set(v.id,{ev,ds:D[i]});if(b.offers.size===2)pending.set(bKey(b),{type:"b",x:b})}
   }else{
    let b=x;if(b.offers.size!==2||!B.has(b.id))throw Error("bad bond");let [a,c]=b.ends.map(v=>b.offers.get(v));let ev="b"+b.n;
    H.set(k,{id:ev,kind:type,g:b.g,parents:[...new Set([b.birth,a.ev,c.ev].filter(Boolean))].sort(),inputs:[a,c].map(o=>o.ds.map(v=>v.n)).sort((a,b)=>cmp(a[0],b[0]))});retiredB++;if(b.roots[0]!==b.roots[1])crossJoins++;B.delete(b.id);
    for(let l=0;l<2;l++)bond(a.ds[l],c.ds[l],b.g+1,ev);
   }
   if(H.size>20000)throw Error("resource guard")
  }
  if([...V.values()].some(v=>v.g!==R||v.bonds.length!==v.need)||[...B.values()].some(b=>b.g!==R||b.offers.size))throw Error("incomplete cut");
  let vol=Object.fromEntries(roots.map(r=>[r,0])),boundary=Object.fromEntries(roots.map(r=>[r,0])),cross={};for(let v of V.values())vol[v.root]++;
  for(let b of B.values()){let [a,c]=b.roots;if(a!==c){boundary[a]++;boundary[c]++;let k=pair(a,c);cross[k]=(cross[k]??0)+1}}
  for(let r of roots)if(boundary[r]!==seed.filter(e=>e.includes(r)).length*2**R)throw Error("boundary");
  if(Object.keys(cross).length!==seed.length||seed.some(([a,b])=>cross[pair(a,b)]!==2**R))throw Error("cross pair");
  if(seed.length===6&&roots.length===4&&roots.some(r=>vol[r]!==6**R))throw Error("K4 volume");
  let q=[V.keys().next().value],seen=new Set(q);for(let i=0;i<q.length;i++)for(let id of V.get(q[i]).bonds){let e=B.get(id),w=e.ends.find(v=>v!==q[i]);if(!seen.has(w)){seen.add(w);q.push(w)}}if(seen.size!==V.size)throw Error("disconnected");
  let vs=[...V.values()].map(v=>[v.n,v.root,v.g,v.birth,v.need,v.bonds.map(id=>B.get(id).n).sort()]).sort((a,b)=>cmp(a[0],b[0]));
  let bs=[...B.values()].map(b=>[b.n,b.g,b.birth]).sort((a,b)=>cmp(a[0],b[0]));
  let hs=[...H.values()].sort((a,b)=>cmp(a.id,b.id));
  return {sig:JSON.stringify([vs,bs,hs]),metrics:{V:V.size,B:B.size,cycles:B.size-V.size+1,retiredV,retiredB,crossJoins,maxD,maxEnabled,adjacent,vol,boundary,cross,connected:true,mincut:flowOn?mincut(V,B,seed[0]):null}};
 }
 let records=[];
 for(let [name,seed] of Object.entries(seeds))for(let R=1;R<=3;R++){let base;for(let [mode,rs] of [["asc",0],["desc",0],["random",17],["random",29]])for(let relabel of [false,true]){
  let r=run(seed,R,mode,rs,relabel?[13,7,19,2]:[0,1,2,3],mode==="asc"&&!relabel);if(base===undefined)base=r.sig;if(base!==r.sig)throw Error("signature mismatch "+name+R+mode+relabel);records.push({name,R,mode,rs,relabel,match:true,...r.metrics});
 }}
 return {records};
}
if (typeof require !== 'undefined' && require.main === module) console.log(JSON.stringify(run31E(),null,2));
