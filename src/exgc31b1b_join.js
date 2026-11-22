function runGC31B1B(){
 const seeds={single:[[0,1]],path:[[0,1],[1,2],[2,3]],triangle:[[0,1],[1,2],[2,0]],cycle:[[0,1],[1,2],[2,3],[3,0]],star:[[0,1],[0,2],[0,3]]};
 const pair=(a,b)=>a<b?`(${a}|${b})`:`(${b}|${a})`;
 const compare=(a,b)=>a<b?-1:a>b?1:0;
 function run(seed,R,mode,randomSeed,perm){
  let edges=new Map, vertices=new Map, incident=new Map, pending=new Map, events=new Map, grade=Array(R+1).fill(0), maxEnabled=0, sharedOverlap=false, initialDegrees;
  let rng=randomSeed>>>0;function rand(n){rng=(Math.imul(1664525,rng)+1013904223)>>>0;return rng%n;}
  function addVertex(v,n){if(vertices.has(v))throw Error("duplicate vertex");vertices.set(v,n);incident.set(v,new Set);}
  function offer(e,f){if(e.id===f.id)return;let k=pair(e.id,f.id),r=1+Math.max(e.rank,f.rank);if(r>R||pending.has(k)||events.has(k))return;let common=e.ends.filter(v=>f.ends.includes(v));if(common.length!==1)return;pending.set(k,{id:k,e,f,r,outer:[e.ends.find(v=>v!==common[0]),f.ends.find(v=>v!==common[0])]});}
  function addEdge(e){if(edges.has(e.id))throw Error("duplicate edge");edges.set(e.id,e);for(let v of e.ends){for(let f of incident.get(v))offer(e,edges.get(f));incident.get(v).add(e.id);} }
  for(let v of [...new Set(seed.flat())])addVertex("s"+perm[v],"s"+v);
  for(let [a,b] of seed){let ends=["s"+perm[a],"s"+perm[b]],id="i"+pair(...ends);addEdge({id,n:"i"+pair("s"+a,"s"+b),ends,rank:0,producer:null,np:null});}
  initialDegrees=[...vertices.keys()].map(v=>[vertices.get(v),incident.get(v).size]).sort();
  while(pending.size){
   maxEnabled=Math.max(maxEnabled,pending.size);
   if(!sharedOverlap){let seen=new Set;for(let p of pending.values()){if(seen.has(p.e.id)||seen.has(p.f.id)){sharedOverlap=true;break}seen.add(p.e.id);seen.add(p.f.id)}}
   let ks=[...pending.keys()].sort(compare), k=mode==="asc"?ks[0]:mode==="desc"?ks[ks.length-1]:ks[rand(ks.length)];
   let p=pending.get(k);pending.delete(k);
   let norm=pair(p.e.n,p.f.n),w="v"+k,nw="v"+norm;
   events.set(k,{n:norm,parents:[p.e.np,p.f.np].filter(x=>x!==null).sort(compare),rank:p.r});grade[p.r]++;
   addVertex(w,nw);
   for(let v of p.outer)addEdge({id:"o"+pair(w,v),n:"o"+pair(nw,vertices.get(v)),ends:[w,v],rank:p.r,producer:k,np:norm});
   if(events.size>200000)throw Error("resource guard");
  }
  let vs=[...vertices.values()].sort(compare);
  let es=[...edges.values()].map(e=>[e.n,e.ends.map(v=>vertices.get(v)).sort(compare),e.rank,e.np]).sort((a,b)=>compare(a[0],b[0]));
  let ev=[...events.values()].sort((a,b)=>compare(a.n,b.n));
  let degrees=[...incident.values()].map(s=>s.size), n=vertices.size,diameter=null,connected=null;
  if(n<=2000){diameter=0;connected=true;for(let v of vertices.keys()){let dist=new Map([[v,0]]),q=[v];for(let i=0;i<q.length;i++){let x=q[i];for(let id of incident.get(x)){let y=edges.get(id).ends.find(y=>y!==x);if(!dist.has(y)){dist.set(y,dist.get(x)+1);q.push(y)}}}if(q.length!==n)connected=false;diameter=Math.max(diameter,...dist.values());}}
  else {let first=vertices.keys().next().value,q=[first],seen=new Set(q);for(let i=0;i<q.length;i++)for(let id of incident.get(q[i])){let y=edges.get(id).ends.find(y=>y!==q[i]);if(!seen.has(y)){seen.add(y);q.push(y)}}connected=seen.size===n;}
  return {signature:JSON.stringify([vs,es,ev]),metrics:{events:events.size,grade,vertices:n,edges:edges.size,maxDegree:Math.max(...degrees),meanDegree:2*edges.size/n,diameter,connected,maxEnabled,sharedOverlap,initialDegrees,seedDegrees:[...vertices.keys()].filter(v=>v.startsWith("s")).map(v=>[vertices.get(v),incident.get(v).size]).sort()}};
 }
 let records=[];
 for(let [name,seed] of Object.entries(seeds))for(let R=1;R<=3;R++){
  let base;for(let [mode,s] of [["asc",0],["desc",0],["random",17],["random",29]])for(let relabel of [false,true]){
   let t=run(seed,R,mode,s,relabel?[13,7,19,2]:[0,1,2,3]);if(!base)base=t.signature;if(t.signature!==base)throw Error("schedule mismatch "+name+R);
   records.push({name,R,mode,s,relabel,match:true,...t.metrics});
  }
 }
 // Destructive variant: path edges 0,1,2,3. Delete middle pair, or end pairs.
 function destructive(order){let live=new Set([0,1,2,3]),events=0;for(let i of order){if(!live.has(i)||!live.has(i+1))throw Error("bad control");live.delete(i);live.delete(i+1);events++;}let terminal=![...live].some(i=>live.has(i+1));return {events,terminal,remaining:[...live]};}
 return {records,negative:{central:destructive([1]),outer:destructive([0,2])}};
}
if (typeof require !== 'undefined' && require.main === module) console.log(JSON.stringify(runGC31B1B(),null,2));
