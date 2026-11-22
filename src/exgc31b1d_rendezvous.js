function run31D(){
 const seeds={single:[[0,1]],path:[[0,1],[1,2],[2,3]],triangle:[[0,1],[1,2],[2,0]],cycle:[[0,1],[1,2],[2,3],[3,0]],star:[[0,1],[0,2],[0,3]]};
 const pair=(a,b)=>a<b?"("+a+"|"+b+")":"("+b+"|"+a+")";
 const cmp=(a,b)=>a<b?-1:a>b?1:0;
 function run(seed,R,arm,mode,rs,perm){
  let V=new Map,B=new Map,P=new Map,H=new Map,retired=0,split=0,joins=0,cross=0,crossAny=0,adjacent=false,peakQ=1,peakD=0,splitRoots={},gradeCounts=Array(R+1).fill(0),rng=rs>>>0;
  let roots=[...new Set(seed.flat())].sort((a,b)=>a-b);
  function addV(id,n,root,birth,rank){if(V.has(id)||H.has(id))throw Error("duplicate vertex");let v={id,n,root,birth,rank,ports:[]};V.set(id,v);return v;}
  function tok(prod,rank,root,bit){return {prod,rank,root,bit}}
  function bond(a,b,prod,rank,root,bit,markA=false){
   let id=pair(a.id,b.id),n=pair(a.n,b.n);if(B.has(id))throw Error("duplicate bond");
   let pa={id:id+"@"+a.id,n:n+"@"+a.n,owner:a.id,q:[tok(prod,rank,root,markA?1:bit)]},pb={id:id+"@"+b.id,n:n+"@"+b.n,owner:b.id,q:[tok(prod,rank,root,bit)]};pa.other=pb.id;pb.other=pa.id;
   P.set(pa.id,pa);P.set(pb.id,pb);B.set(id,{n,ends:[pa.id,pb.id]});a.ports.push(pa.id);b.ports.push(pb.id);
  }
  for(let r of roots)addV("s"+perm[r],"s"+r,r,null,0);
  for(let i=0;i<seed.length;i++){let [a,b]=seed[i];bond(V.get("s"+perm[a]),V.get("s"+perm[b]),null,0,null,0,arm==="R1"&&i===0)}
  function ready(v){let ps=v.ports.map(p=>P.get(p));if(ps.some(p=>!p.q.length))return null;let rank=1+Math.max(v.rank,...ps.map(p=>p.q[0].rank));return rank<=R?{v,ps,rank}:null}
  while(true){
   let list=[...V.values()].map(ready).filter(Boolean).sort((a,b)=>cmp(a.v.id,b.v.id));if(!list.length)break;
   if(!adjacent){let set=new Set(list.map(x=>x.v.id));adjacent=list.some(x=>x.ps.some(p=>set.has(P.get(p.other).owner)))}
   let x=mode==="asc"?list[0]:mode==="desc"?list.at(-1):list[(rng=(Math.imul(1664525,rng)+1013904223)>>>0)%list.length];
   let {v,ps,rank}=x;let inputs=ps.map(p=>({port:p.n,...p.q.shift()})).sort((a,b)=>cmp(a.port,b.port));let bit=inputs.some(t=>t.bit)?1:0;
   let parents=[...new Set([v.birth,...inputs.map(t=>t.prod)].filter(t=>t!==null))].sort();
   let tr=new Set(inputs.filter(t=>t.prod!==null).map(t=>t.root));if(parents.length>=2)joins++;if(tr.size>=2)cross++;if(inputs.some(t=>t.prod!==null&&t.root!==v.root))crossAny++;
   let doSplit=arm==="S"||bit===1;let d=ps.length;
   H.set(v.id,{n:v.n,root:v.root,rank,birth:v.birth,inputs,parents,bit,doSplit,d});gradeCounts[rank]++;retired++;V.delete(v.id);
   if(doSplit&&d>1){split++;splitRoots[v.root]=(splitRoots[v.root]??0)+1}
   let daughters=[];
   if(doSplit){for(let p of ps){let w=addV(v.id+"<"+p.id+">",v.n+"<"+p.n+">",v.root,v.n,rank);w.ports.push(p.id);p.owner=w.id;daughters.push(w)}}
   else {let w=addV(v.id+"c",v.n+"c",v.root,v.n,rank);for(let p of ps){w.ports.push(p.id);p.owner=w.id}daughters.push(w)}
   for(let p of ps){let opp=P.get(p.other);opp.q.push(tok(v.n,rank,v.root,bit));peakQ=Math.max(peakQ,opp.q.length);if(opp.q.length>2)throw Error("queue bound")}
   for(let i=0;i<daughters.length;i++)for(let j=i+1;j<daughters.length;j++)bond(daughters[i],daughters[j],v.n,rank,v.root,bit);
   for(let w of daughters){peakD=Math.max(peakD,w.ports.length);if(w.ports.length>3||w.ports.length<1)throw Error("degree bound")}
   if(H.size>5000)throw Error("guard");
  }
  let ev=[...H.values()].sort((a,b)=>cmp(a.n,b.n));let vs=[...V.values()].map(v=>[v.n,v.root,v.birth,v.rank,v.ports.map(p=>P.get(p).n).sort()]).sort((a,b)=>cmp(a[0],b[0]));
  let bs=[...B.values()].map(b=>[b.n,b.ends.map(p=>{let e=P.get(p);return [e.n,V.get(e.owner).n,e.q]}).sort((a,b)=>cmp(a[0],b[0]))]).sort((a,b)=>cmp(a[0],b[0]));
  let adj=new Map([...V.keys()].map(v=>[v,[]]));for(let b of B.values()){let [a,c]=b.ends.map(p=>P.get(p).owner);if(a===c)throw Error("loop");adj.get(a).push(c);adj.get(c).push(a)}
  let q=[V.keys().next().value],seen=new Set(q);for(let i=0;i<q.length;i++)for(let w of adj.get(q[i]))if(!seen.has(w)){seen.add(w);q.push(w)}if(seen.size!==V.size)throw Error("disconnect");
  let expectedV=roots.length,expectedB=seed.length;for(let e of ev)if(e.doSplit){expectedV+=e.d-1;expectedB+=e.d*(e.d-1)/2}
  if(V.size!==expectedV||B.size!==expectedB||retired!==H.size)throw Error("count");
  if(arm==="R0"&&(split||V.size!==roots.length||B.size!==seed.length||ev.some(e=>e.bit)))throw Error("null");
  let finalD=Math.max(...[...V.values()].map(v=>v.ports.length));peakD=Math.max(peakD,finalD,...roots.map(r=>seed.filter(e=>e.includes(r)).length));
  return {sig:JSON.stringify([ev,vs,bs]),metrics:{events:H.size,retired,gradeCounts,V:V.size,B:B.size,cycles:B.size-V.size+1,peakD,peakQ,split,splitRoots,joins,cross,crossAny,adjacent,connected:true,seedBondsRetained:seed.every(([a,b])=>B.has(pair("s"+perm[a],"s"+perm[b]))) }};
 }
 let records=[];
 for(let arm of ["S","R0","R1"])for(let [name,seed] of Object.entries(seeds))for(let R=1;R<=4;R++){let base;for(let [mode,rs] of [["asc",0],["desc",0],["random",17],["random",29]])for(let relabel of [false,true]){
   let r=run(seed,R,arm,mode,rs,relabel?[13,7,19,2]:[0,1,2,3]);if(base===undefined)base=r.sig;if(base!==r.sig)throw Error("signature mismatch "+arm+name+R+mode+relabel);records.push({arm,name,R,mode,rs,relabel,match:true,...r.metrics});
 }}
 return {records};
}
if (typeof require !== 'undefined' && require.main === module) console.log(JSON.stringify(run31D(),null,2));
