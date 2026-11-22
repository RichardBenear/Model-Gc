function run31F(){
 const configs=[["ring4",4,1,-1],["ring5",5,1,-1],["ring6",6,1,-1],["chord7",7,1,2],["two3",6,1,-1]];
 const cmp=(a,b)=>a<b?-1:a>b?1:0;
 const pair=(a,b)=>a<b?a+"|"+b:b+"|"+a;
 function flow(vs,bs){
  let ix=new Map(vs.map((v,i)=>[v.id,i])),N=vs.length,S=N,T=N+1,adj=Array.from({length:N+2},()=>[]),outside=new Set;
  function edge(a,b,c){if(a===b)return;let e={b,c,rev:adj[b].length},f={b:a,c:0,rev:adj[a].length};adj[a].push(e);adj[b].push(f)}
  for(let b of bs){let a=ix.get(b.a.id),c=ix.get(b.b.id);edge(a,c,1);edge(c,a,1);if(b.a.root!==b.b.root){outside.add(a);outside.add(c)}}
  let L=[],R=[];for(let i=0;i<N;i++)if(!outside.has(i)){if(vs[i].root===0)L.push(i);if(vs[i].root===1)R.push(i)}
  if(!L.length||!R.length)return {value:null,left:L.length,right:R.length};
  for(let a of L)edge(S,a,bs.length+1);for(let b of R)edge(b,T,bs.length+1);
  let value=0;while(true){let prev=Array(N+2).fill(null),q=[S];prev[S]=[-1,-1];for(let j=0;j<q.length&&!prev[T];j++){let a=q[j];for(let i=0;i<adj[a].length;i++){let e=adj[a][i];if(e.c>0&&!prev[e.b]){prev[e.b]=[a,i];q.push(e.b)}}}if(!prev[T])break;let d=Infinity;for(let b=T;b!==S;){let[a,i]=prev[b];d=Math.min(d,adj[a][i].c);b=a}for(let b=T;b!==S;){let[a,i]=prev[b],e=adj[a][i];e.c-=d;adj[b][e.rev].c+=d;b=a}value+=d}
  return {value,left:L.length,right:R.length};
 }
 function run(conf,R,arm,mode,rs,relabel,flowOn){
  let [name,N,step0,step1]=conf,W=arm==="P"?["0","1","00","01","10","11"]:["0","1"],V=new Map,B=[],H=[],pending=new Map,requests=new Set,rng=rs>>>0,joins=0,relay=0,split=0,interRoot=0,overlap=false,peak=0,retired=0;
  function vertex(id,n,root,g,birth){let v={id,n,root,g,birth,out:[null,null],in:[null,null],done:false,waiting:[],children:new Map};if(V.has(id))throw Error("duplicate vertex");V.set(id,v);return v}
  function ready(v){if(!v.done&&v.g<R&&v.out.every(Boolean)&&v.in.every(Boolean))pending.set("v"+v.id,{type:"v",v})}
  function bond(a,b,c,prod){if(a.out[c]||b.in[c]||a.g!==b.g)throw Error("port conflict");let x={a,b,c,prod};a.out[c]=x;b.in[c]=x;B.push(x);ready(a);ready(b);}
  function parentList(xs){return [...new Set(xs.filter(Boolean))].sort()}
  let seed=[];for(let r=0;r<N;r++)seed.push(vertex("s"+(relabel?(N-1-r)*3+7:r),"s"+r,r,0,null));
  function target(r,step){if(name==="two3"){let base=r-r%3;return base+(r%3+step+3)%3}return (r+step+N)%N}
  let initial=new Set;for(let r=0;r<N;r++)for(let [c,s] of [[0,step0],[1,step1]]){let t=target(r,s);bond(seed[r],seed[t],c,null);if(r!==t)initial.add(pair(r,t))}
  function submit(req,v,type){let key=(type==="relay"?"r":"j")+req.key;if(v.done)pending.set(key,{type,req,v});else v.waiting.push({key,type,req,v});}
  while(pending.size){
   peak=Math.max(peak,pending.size);if(!overlap){let set=new Set([...pending.values()].filter(e=>e.type==="v").map(e=>e.v.id));overlap=[...pending.values()].some(e=>e.type==="v"&&e.v.out.some(b=>set.has(b.b.id)&&b.b!==e.v))}
   let ks=[...pending.keys()].sort(cmp),k=mode==="asc"?ks[0]:mode==="desc"?ks.at(-1):ks[(rng=(Math.imul(1664525,rng)+1013904223)>>>0)%ks.length],e=pending.get(k);pending.delete(k);
   let v=e.v;
   if(e.type==="v"){
    if(v.done||!v.out.every(Boolean)||!v.in.every(Boolean))throw Error("bad ready");
    v.done=true;retired++;split++;let ev="v"+v.n;
    H.push({id:ev,g:v.g,type:"v",parents:parentList([v.birth,...v.out.map(b=>b.prod),...v.in.map(b=>b.prod)])});
    for(let w of W)v.children.set(w,vertex(v.id+"/"+w,v.n+"/"+w,v.root,v.g+1,ev));
    for(let i=0;i<W.length;i++)bond(v.children.get(W[i]),v.children.get(W[(i+1)%W.length]),1,ev);
    for(let p of v.waiting)pending.set(p.key,p);v.waiting=[];
    for(let w of W){let x=v.out[+w[0]],req={src:v,w,key:v.id+"/"+w,n:v.n+"/"+w,prev:ev,read:x.prod};if(requests.has(req.key))throw Error("duplicate request");requests.add(req.key);submit(req,x.b,w.length===2?"relay":"join")}
   }else if(e.type==="relay"){
    let req=e.req;if(!v.done)throw Error("relay not ready");let ev="r"+req.n;let x=v.out[+req.w[1]];
    H.push({id:ev,g:req.src.g,type:"relay",parents:parentList([req.prev,"v"+v.n,req.read,x.prod])});relay++;req.prev=ev;req.read=x.prod;submit(req,x.b,"join");
   }else{
    let req=e.req;if(!v.done||!requests.delete(req.key))throw Error("join not ready");let ev="j"+req.n;H.push({id:ev,g:req.src.g,type:"join",parents:parentList([req.prev,"v"+v.n,req.read])});joins++;if(req.src.root!==v.root)interRoot++;bond(req.src.children.get(req.w),v.children.get(req.w),0,ev);
   }
   if(H.length>30000)throw Error("cap")
  }
  let vs=[...V.values()].filter(v=>!v.done),bs=B.filter(b=>!b.a.done&&!b.b.done);
  if(requests.size||vs.some(v=>v.g!==R||!v.out.every(Boolean)||!v.in.every(Boolean))||vs.length!==N*W.length**R||bs.length!==2*vs.length)throw Error("incomplete cut");
  if([...V.values()].some(v=>v.waiting.length))throw Error("unprocessed wait");
  let rootPairs=new Set,boundary=Object.fromEntries(seed.map(v=>[v.root,0])),vol={...boundary},loops=0;for(let v of vs)vol[v.root]++;
  for(let b of bs){if(b.a===b.b)loops++;if(b.a.root!==b.b.root){rootPairs.add(pair(b.a.root,b.b.root));boundary[b.a.root]++;boundary[b.b.root]++}}
  let newPairs=[...rootPairs].filter(p=>!initial.has(p)).sort(),missing=[...initial].filter(p=>!rootPairs.has(p));
  if(missing.length||arm==="S"&&newPairs.length||name==="two3"&&[...rootPairs].some(p=>{let[a,b]=p.split("|").map(Number);return Math.floor(a/3)!==Math.floor(b/3)}))throw Error("root invariant");
  let seen=new Set,components=0;for(let v of vs)if(!seen.has(v.id)){components++;let q=[v];seen.add(v.id);for(let i=0;i<q.length;i++)for(let b of [...q[i].out,...q[i].in]){let w=b.a===q[i]?b.b:b.a;if(!seen.has(w.id)){seen.add(w.id);q.push(w)}}}
  if(components!==(name==="two3"?2:1))throw Error("components");
  let vv=vs.map(v=>[v.n,v.root,v.g,v.birth]).sort((a,b)=>cmp(a[0],b[0]));
  let bb=bs.map(b=>[b.a.n,b.b.n,b.c,b.prod]).sort((a,b)=>cmp(a[0]+":"+a[2],b[0]+":"+b[2]));H.sort((a,b)=>cmp(a.id,b.id));
  return {sig:JSON.stringify([vv,bb,H]),metrics:{V:vs.length,B:bs.length,split,relay,joins,interRoot,retired,overlap,peak,maxDegree:4,loops,components,rootPairs:rootPairs.size,initialPairs:initial.size,newPairs,boundary,vol,mincut:flowOn?flow(vs,bs):null}};
 }
 let records=[];
 for(let arm of ["P","S"])for(let conf of configs)for(let R=1;R<=3;R++){let base;for(let [mode,rs] of [["asc",0],["desc",0],["random",17],["random",29]])for(let relabel of [false,true]){let r=run(conf,R,arm,mode,rs,relabel,mode==="asc"&&!relabel);if(base===undefined)base=r.sig;if(base!==r.sig)throw Error("signature mismatch "+arm+conf[0]+R);records.push({arm,name:conf[0],R,mode,rs,relabel,match:true,...r.metrics})}}
 return {records}
}
if (typeof require !== 'undefined' && require.main === module) console.log(JSON.stringify(run31F(),null,2));
