interval = 8

f = open('cube.raw', 'w')

for i in xrange(0,64):
    for j in xrange(0,64):
        for k in xrange(0,64):
            if i > 16 and i < 48 and \
               j > 16 and j < 48 and \
               k > 16 and j < 48:
                   if (not ((k / interval) % 2)) == (j / interval) % 2:
                       f.write("%c" % 255)
                   else:
                       f.write("%c" % 254)
            else:
                f.write("%c" % 253)

    print("%d\n" % i)
