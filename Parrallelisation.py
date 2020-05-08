mblocklistnum = []
for i in range(len(mblocklist)):
    mblocklistnum.append([i, mblocklist[i]])
p = Pool(multiprocessing.cpu_count())
results = p.map(partial(pool_write_microalignment,targetdata=targetdata,extendedsourcedata=extendedsourcedata,nbinitialsource=nbinitialsource,all_ids=all_ids,msamethod=msamethod), mblocklistnum)

def pool_write_microalignment(mblocknum,targetdata,extendedsourcedata,nbinitialsource,all_ids,msamethod)
#dernier param = que'est-ce qui change
