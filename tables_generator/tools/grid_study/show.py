import json,sys
for f in sys.argv[1:]:
  for l in open(f):
    d=json.loads(l)
    s=' '.join('%g:%s[Zh%.3f Zdr%.3f r%.4f K%.3f A%.3f]'%(float(k),'OK' if v['ok'] else '--',v['Zh'],v['Zdr'],v['rhohv'],v['Kdp'],v['Ah']) for k,v in d['steps'].items())
    print('%-2s %-2s cur %-6g ref %-6g| %s'%(d['band'],d['type'],d['current'],d['ref'],s))
