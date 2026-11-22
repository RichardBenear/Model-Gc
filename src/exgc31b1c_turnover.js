function run31C(){
 const seeds={single:[[0,1]],path:[[0,1],[1,2],[2,3]],triangle:[[0,1],[1,2],[2,0]],cycle:[[0,1],[1,2],[2,3],[3,0]],star:[[0,1],[0,2],[0,3]]};
 const pair=(a,b)=>a<b?"("+a+"|"+b+")":"("+b+"|"+a+")";
 function run(seed,R,arm,mode,rs,perm){
  let verts=new Map,adj=new Map,live=new Map,pending=new Map,history=new Map,retired=new Set,grades=Array(R+1).fill(0),maxEnabled=0,overlap=false;
  let rng=rs>>>0;const rand=n=>{rng=(Math.imul(1664525,rng)+1013904223)>>>0;return rng%n};
  function vertex(v,n){if(verts.has(v))throw Error("duplicate vertex");verts.set(v,n);adj.set(v,new Set);}
  function edge(a,b,rank,producer,np){let id=pair(a,b),n=pair(verts.get(a),verts.get(b));if(live.has(id)||retired.has(id))throw Error("duplicate instance");let e={id,n,a,b,rank,producer,np};live.set(id,e);adj.get(a).add(b);adj.get(b).add(a);if(rank<R)pending.set(id,e);}
  let roots=[...new Set(seed.flat())].sort((a,b)=>a-b);
  for(let v of roots)vertex("s"+perm[v],"s"+v);
  for(let [a,b] of seed)edge("s"+perm[a],"s"+perm[b],0,null,null);
  let initialMax=Math.max(...[...adj.values()].map(x=>x.size)),bound=Math.max(initialMax,arm==="S"?2:3),peak=initialMax;
  function degree(v){let d=adj.get(v).size;peak=Math.max(peak,d);if(d>bound)throw Error("degree bound");}
  while(pending.size){
   maxEnabled=Math.max(maxEnabled,pending.size);if(!overlap){let seen=new Set;for(let e of pending.values()){if(seen.has(e.a)||seen.has(e.b)){overlap=true;break}seen.add(e.a);seen.add(e.b)}}
   let keys=[...pending.keys()].sort(),id=mode==="asc"?keys[0]:mode==="desc"?keys.at(-1):keys[rand(keys.length)],e=pending.get(id);pending.delete(id);live.delete(id);retired.add(id);adj.get(e.a).delete(e.b);adj.get(e.b).delete(e.a);
   let [a,b]=[e.a,e.b].sort(),rank=e.rank+1,k=e.id,n=e.n,added=[];
   history.set(k,{n,rank,parent:e.np});grades[rank]++;
   function fresh(tag,ntag){let v="v"+k+tag;vertex(v,"v"+n+ntag);added.push(v);return v}
   function add(x,y){edge(x,y,rank,k,n)}
   if(arm==="S"){let w=fresh("m","m");add(a,w);add(w,b);}
   else {let p=fresh("p","t"+verts.get(a)),q=fresh("q","t"+verts.get(b)),r=fresh("r","r"),s=fresh("s","s");add(a,p);add(q,b);add(p,r);add(r,q);add(p,s);add(s,q);}
   for(let v of [a,b,...added])degree(v);
  }
  let dist=[];
  for(let i=0;i<roots.length;i++){let a="s"+perm[roots[i]],q=[a],seen=new Map([[a,0]]);for(let j=0;j<q.length;j++){for(let v of adj.get(q[j]))if(!seen.has(v)){seen.set(v,seen.get(q[j])+1);q.push(v)}}if(q.length!==verts.size)throw Error("disconnected");for(let j=i+1;j<roots.length;j++)dist.push([roots[i],roots[j],seen.get("s"+perm[roots[j]])]);}
  // Compute independent seed distances directly.
  let sd=new Map(roots.map(x=>[x,[]]));for(let [a,b] of seed){sd.get(a).push(b);sd.get(b).push(a)}
  let L=arm==="S"?2:4,K=arm==="S"?2:6,Q=arm==="S"?1:4;
  for(let [a,b,d] of dist){let queue=[a],D=new Map([[a,0]]);for(let i=0;i<queue.length;i++)for(let v of sd.get(queue[i]))if(!D.has(v)){D.set(v,D.get(queue[i])+1);queue.push(v)}if(d!==D.get(b)*L**R)throw Error("distance scale");}
  let m=seed.length,n=roots.length,expected=m*(K**R-1)/(K-1);
  if(history.size!==expected||live.size!==m*K**R||verts.size!==n+Q*expected||retired.size!==expected||[...live.values()].some(x=>x.rank!==R))throw Error("count/turnover");
  let vs=[...verts.values()].sort(),es=[...live.values()].map(x=>[x.n,x.rank,x.np]).sort((a,b)=>a[0]<b[0]?-1:a[0]>b[0]?1:0),ev=[...history.values()].sort((a,b)=>a.n<b.n?-1:a.n>b.n?1:0);
  let sum=JSON.stringify([vs,es,ev]);return {sum,metrics:{events:history.size,retired:retired.size,grades,V:verts.size,E:live.size,cycleRank:live.size-verts.size+1,peak,finalMax:Math.max(...[...adj.values()].map(x=>x.size)),bound,maxEnabled,overlap,dist}};
 }
 let records=[];
 for(let arm of ["S","D"])for(let [name,seed] of Object.entries(seeds))for(let R=1;R<=3;R++){let base;for(let [mode,rs] of [["asc",0],["desc",0],["random",17],["random",29]])for(let relabel of [false,true]){
  let r=run(seed,R,arm,mode,rs,relabel?[13,7,19,2]:[0,1,2,3]);if(base===undefined)base=r.sum;if(base!==r.sum)throw Error("signature mismatch "+arm+name+R);records.push({arm,name,R,mode,rs,relabel,match:true,...r.metrics});
 }}
 return {records};
}
if (typeof require !== 'undefined' && require.main === module) console.log(JSON.stringify(run31C(),null,2));
