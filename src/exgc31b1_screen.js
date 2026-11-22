function screen(csv) {
 const rel=(a,b,c)=>[a,Math.min(b,c),Math.max(b,c)].join(",");
 const lhs=([x,y,z,u,v])=>[rel(x,y,z),rel(y,z,x),rel(x,u,v),rel(v,y,z)];
 const key=a=>a.join(",");
 function gates(s) {
  let by=new Map;for(let k of s.keys()){let [a,b,c]=k.split(",").map(Number);if(!by.has(a))by.set(a,[]);by.get(a).push([b,c]);}
  let out=[];
  for(let [x,ps] of by)for(let [b,c] of ps)for(let [y,z] of [[b,c],[c,b]])if(s.has(rel(y,z,x)))for(let [d,e] of ps)for(let [u,v] of [[d,e],[e,d]])if(new Set([x,y,z,u,v]).size===5&&s.has(rel(v,y,z)))out.push([x,y,z,u,v]);
  return out.sort((a,b)=>key(a).localeCompare(key(b)));
 }
 function fire(s,g,output,w,event) {
  let req=lhs(g),t=new Map(s),dep=0;
  for(let i=0;i<4;i++)if(s.get(req[i])===1)dep|=i===3?2:1;
  for(let k of req.slice(0,3))t.delete(k);
  let vals=[0,...g.slice(1),w];
  for(let [a,b,c] of output){let k=rel(vals[a],vals[b],vals[c]);if(!t.has(k))t.set(k,event);}
  return {s:t,dep};
 }
 function sig(s,w,swap=false){return [...s.keys()].map(k=>{let q=k.split(",").map(Number).map(x=>swap?(x===w?w+1:x===w+1?w:x):x);return rel(...q)}).sort().join(";")}
 let contexts=[],base=[1,2,3,4,5],a=lhs(base);
 function gen(b,used,next){
  if(b.length<5){for(let k=1;k<=5;k++)if(!used.has(k)){let n=new Set(used);n.add(k);gen([...b,k],n,next)}gen([...b,next],used,next+1);return}
  if(!used.size||key(b)===key(base))return;
  let q=lhs(b);if(!a.slice(0,3).some(k=>q.includes(k))&&!q.slice(0,3).some(k=>a.includes(k)))return;
  let s=new Map([...a,...q].map(k=>[k,0]));contexts.push({b,s,w:next});
 }
 gen([],new Set,6);
 let rows=csv.trim().split(/\r?\n/).slice(1).map(l=>{let [id]=l.split(",");let output=[...l.matchAll(/\((\d+), (\d+), (\d+)\)/g)].map(m=>m.slice(1).map(Number));return {id:+id,output}});
 if(rows.length!==80||new Set(rows.map(r=>r.id)).size!==80||rows.some(r=>r.output.length!==8))throw Error("bad population");
 let results=[];
 for(let {id,output} of rows){
  let disjoint=new Map([...a,...lhs([6,7,8,9,10])].map(k=>[k,0]));
  let d1=fire(fire(disjoint,base,output,11,1).s,[6,7,8,9,10],output,12,2);
  let d2=fire(fire(disjoint,[6,7,8,9,10],output,11,1).s,base,output,12,2);
  if(d1.dep||d2.dep||sig(d1.s,11)!==sig(d2.s,11,true))throw Error("disjoint control failed "+id);
  let local=0,terminalBoth=0,diamond=0,first=null,strong=null,tested=0;
  for(let {b,s,w} of contexts){
   tested++;
   let A=fire(s,base,output,w,1),B=fire(s,b,output,w,1);let ga=gates(A.s),gb=gates(B.s);
   if((ga.length===0)!==(gb.length===0)){strong={b,w,seed:[...s.keys()].sort(),A:[...A.s.keys()].sort(),B:[...B.s.keys()].sort(),nextA:ga,nextB:gb};break}
   if(!ga.length&&!gb.length){terminalBoth++;if(sig(A.s,w)!==sig(B.s,w)){local++;if(!first)first={b,w,kind:"both-terminal-rooted-state-difference"}}continue}
   let set=new Set;for(let g of ga){let t=fire(A.s,g,output,w+1,2);set.add(t.dep+"|"+sig(t.s,w));set.add(t.dep+"|"+sig(t.s,w,true))}
   let ok=gb.some(g=>{let t=fire(B.s,g,output,w+1,2);return set.has(t.dep+"|"+sig(t.s,w))});
   if(ok)diamond++;else{local++;if(!first)first={b,w,kind:"no-two-event-rooted-diamond",nextA:ga.length,nextB:gb.length}}
  }
  results.push({id,tested,local,terminalBoth,diamond,first,strong,disjoint:"pass"})
 }
 return {contexts:contexts.length,results};
}
// Node usage: node src/exgc31b1_screen.js W21_CLOSED_PARTICLE_HITS.csv
if (typeof require !== 'undefined' && require.main === module) { const fs=require('node:fs'); console.log(JSON.stringify(screen(fs.readFileSync(process.argv[2], 'utf8')),null,2)); }
