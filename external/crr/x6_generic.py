import random, collections, itertools, json, copy

def rel(a,b,c):
    if b>c:b,c=c,b
    return (a,b,c)

REL_UNIVERSE=[]
for a in range(5):
    for b,c in itertools.combinations([i for i in range(5) if i!=a],2):
        REL_UNIVERSE.append((a,b,c))

def decode_mask(m):
    return {REL_UNIVERSE[i] for i in range(30) if (int(m)>>i)&1}

class RuleSpec:
    def __init__(self,template_index,template,state_masks,period=None):
        self.template_index=int(template_index)
        self.template=[tuple(map(int,r)) for r in template]
        self.phase_rels=[decode_mask(m) for m in state_masks]
        self.period=period
    def outputs(self,y,z,u,v,w):
        roles={1:y,2:z,3:u,4:v,5:w}
        return [rel(roles[a],roles[b],roles[c]) for a,b,c in self.template]

class Sim:
    def __init__(self,spec,seed=0):
        self.spec=spec; self.rng=random.Random(seed)
        self.rels=set(); self.by_src=collections.defaultdict(set); self.target_sources=collections.defaultdict(set); self.n=5
        for r in [rel(0,1,2),rel(1,2,0),rel(0,3,4),rel(4,1,2)]:self._add(r)
        self.gates_by_x={};self.gates_set=set();self._recompute_xs(set(self.by_src))
    def clone(self):
        s=Sim.__new__(Sim);s.spec=self.spec;s.rng=random.Random();s.rng.setstate(self.rng.getstate());s.rels=set(self.rels)
        s.by_src=collections.defaultdict(set,{k:set(v) for k,v in self.by_src.items()});s.target_sources=collections.defaultdict(set,{k:set(v) for k,v in self.target_sources.items()});s.n=self.n
        s.gates_by_x={k:set(v) for k,v in self.gates_by_x.items()};s.gates_set=set(self.gates_set);return s
    def _add(self,r):
        if r in self.rels:return False
        self.rels.add(r);self.by_src[r[0]].add(r);self.target_sources[r[1]].add(r[0]);self.target_sources[r[2]].add(r[0]);return True
    def _remove(self,r):
        if r not in self.rels:return False
        self.rels.remove(r);self.by_src[r[0]].discard(r)
        for t in (r[1],r[2]):
            if not any(rr[1]==t or rr[2]==t for rr in self.by_src[r[0]]):self.target_sources[t].discard(r[0])
        return True
    def _gates_for_x(self,x):
        out=set();R=self.rels;xrels=self.by_src.get(x,())
        if len(xrels)<2:return out
        for r1 in xrels:
            _,p,q=r1
            for y,z in ((p,q),(q,p)):
                if len({x,y,z})<3:continue
                if rel(y,z,x) not in R:continue
                for r3 in xrels:
                    if r3==r1:continue
                    _,u,v=r3
                    if len({x,y,z,u,v})<5:continue
                    if rel(v,y,z) in R:out.add((x,y,z,u,v))
        return out
    def _affected_xs(self,changed):
        xs=set(changed)
        for s in list(changed):
            for r in self.by_src.get(s,()):xs.update((r[1],r[2]))
            xs.update(self.target_sources.get(s,()))
        return xs
    def _recompute_xs(self,xs):
        for x in xs:
            old=self.gates_by_x.get(x,set());self.gates_set.difference_update(old)
            new=self._gates_for_x(x)
            if new:self.gates_by_x[x]=new;self.gates_set.update(new)
            else:self.gates_by_x.pop(x,None)
    def gate_valid(self,g):
        x,y,z,u,v=g
        return len(set(g))==5 and u<v and rel(x,y,z) in self.rels and rel(y,z,x) in self.rels and rel(x,u,v) in self.rels and rel(v,y,z) in self.rels
    def fire(self,g):
        assert self.gate_valid(g),g
        x,y,z,u,v=g;w=self.n;self.n+=1;changed={x,y,z,u,v,w};xs0=self._affected_xs(changed)
        for r in [rel(x,y,z),rel(y,z,x),rel(x,u,v)]:self._remove(r)
        for r in self.spec.outputs(y,z,u,v,w):self._add(r)
        xs1=self._affected_xs(changed);self._recompute_xs(xs0|xs1);return w
    def step(self):
        if not self.gates_set:return None
        g=self.rng.choice(tuple(self.gates_set));self.fire(g);return g
    def run(self,n):
        k=0
        for _ in range(n):
            if self.step() is None:break
            k+=1
        return k
    def induced(self,S):
        S=set(S);return {r for r in self.rels if set(r)<=S}
    def normalized_rels(self,S):
        ss=sorted(S)
        if len(ss)!=5:return None
        mp={v:i for i,v in enumerate(ss)}
        return {(mp[a],min(mp[b],mp[c]),max(mp[b],mp[c])) for a,b,c in self.induced(ss)}
    def phase(self,S):
        rr=self.normalized_rels(S)
        if rr is None:return None
        for i,p in enumerate(self.spec.phase_rels):
            if rr==p:return i
        return None
    def particle_supports(self,phase=None):
        sups=set()
        for g in self.gates_set:
            S=tuple(sorted(g));p=self.phase(S)
            if p is not None and (phase is None or p==phase):sups.add(S)
        return sups
    def internal_gates(self,S):
        ss=set(S);return [g for g in self.gates_set if set(g)==ss]
    def twosection(self):
        A=collections.defaultdict(set)
        for a,b,c in self.rels:
            A[a].update((b,c));A[b].update((a,c));A[c].update((a,b))
        return A

def successor_support(g,w):
    return tuple(sorted((g[1],g[2],g[3],g[4],w)))

def internal_transition(sim,S,g):
    c=sim.clone();w=c.fire(g);Sp=successor_support(g,w);return c,Sp,c.phase(Sp)

def bfs_multi(A,S):
    S=set(S);d={v:0 for v in S};q=collections.deque(S)
    while q:
        v=q.popleft()
        for u in A.get(v,()):
            if u not in d:d[u]=d[v]+1;q.append(u)
    return d

def setmind(A,S,T):
    d=bfs_multi(A,S);vals=[d.get(t,10**9) for t in T];return min(vals) if vals else 10**9

def load_specs(json_path):
    data=json.load(open(json_path))
    return {int(x['template_index']):RuleSpec(x['template_index'],x['template'],x['state_masks'],x.get('period')) for x in data}

def gates_from_rels(R):
    by=collections.defaultdict(set)
    for r in R:by[r[0]].add(r)
    out=set()
    for x,xrels in by.items():
        if len(xrels)<2:continue
        for r1 in xrels:
            _,p,q=r1
            for y,z in ((p,q),(q,p)):
                if len({x,y,z})<3 or rel(y,z,x) not in R:continue
                for r3 in xrels:
                    if r3==r1:continue
                    _,u,v=r3
                    if len({x,y,z,u,v})<5:continue
                    if rel(v,y,z) in R:out.add((x,y,z,u,v))
    return sorted(out)

def normalize_subset(R,S):
    ss=sorted(S);mp={v:i for i,v in enumerate(ss)}
    return {(mp[a],min(mp[b],mp[c]),max(mp[b],mp[c])) for a,b,c in R if set((a,b,c))<=set(S)}

def abstract_phase_graph(spec):
    out={}
    for i,R0 in enumerate(spec.phase_rels):
        rows=[]
        for g in gates_from_rels(R0):
            x,y,z,u,v=g;w=5;R=set(R0)
            for r in [rel(x,y,z),rel(y,z,x),rel(x,u,v)]:R.discard(r)
            for r in spec.outputs(y,z,u,v,w):R.add(r)
            S2=tuple(sorted((y,z,u,v,w)));rr=normalize_subset(R,S2)
            j=None
            for q,p in enumerate(spec.phase_rels):
                if rr==p:j=q;break
            rows.append((g,j,len(rr)))
        out[i]=rows
    return out
