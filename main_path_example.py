import main_path as main_path
def run():
    table=main_path.impcsv("net.csv")
    G=main_path.generateGraph(table)
    coms=main_path.getCluster(G)
    return_path,filterNode,KnodeList = main_path.process(G,coms)
    main_path.createGephiNodeCSV(coms,return_path,KnodeList,filterNode,G)

if __name__ == "__main__":
    run()
