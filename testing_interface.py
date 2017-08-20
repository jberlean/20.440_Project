
import dbm

c = dbm.open('cache','c')
c['a'] = 'c'
c.close()


d = dbm.open('cache','c')
print d['a']
d.close()
