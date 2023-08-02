import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import itertools

with pl.StringCache():
    methodnames = {"999":"orig+", "99":"orig", "0":"opt", "7":"rk4", "4":"sp2", "5":"sp4", "6":"sp6"}
    datanames = ["test47"]
    datanames2 = ["n="+t[4:] for t in datanames]
    pl.Series(methodnames.values()).cast(pl.Categorical)
    pl.Series(datanames2).cast(pl.Categorical)
    dfer = pl.DataFrame()
    for t in datanames:
        d = pl.read_csv(
                source = "../"+t+"-plus/test-results/rock.txt",
                has_header = False,
                columns = [1,2,3,4,5,6],
                separator = " ",
                dtypes = [pl.Categorical, pl.Float64, pl.Float64, pl.Float64, pl.Float64, pl.Float64])
        d.columns = ["method", "dz", "ih", "ik", "angle", "intensity"]
        d = d.with_columns(pl.lit("n="+t[4:]).alias("data")).with_columns(pl.col("method").apply(lambda x: methodnames[x]).cast(pl.Categorical))
        dfer = dfer.vstack(d)

    #selmethods = ["orig+", "sp4", "sp6", "rk4"]
    selmethods = ["sp4", "sp6", "rk4", "orig+"]
    dfall = pl.DataFrame()
    for s in selmethods:
        dfall = dfall.vstack(dfer.filter(pl.col("method") == s))

    pltcmap = lambda x: "#{:02x}{:02x}{:02x}".format(*tuple(((np.array(plt.get_cmap("Set1").colors)*255)[x].astype(int))))
    cmap = {
            #"orig":(pltcmap(0),"o", "solid"),
            #"opt":(pltcmap(1),".", "solid"),
            "rk4":(pltcmap(2),"*", "solid"),
            "sp4":(pltcmap(3),"+", "solid"),
            "sp6":(pltcmap(4),"x", "solid"),
            #"orig+":(pltcmap(1),"o", "dotted"),
            "orig+":("black","", "solid"),
            }

    dzs = [0.68, 0.33, 0.15]
    ihiks=[(0,0), (1,-1), (-1,1)]#, (2,-2), (-2,2), (12,-12), (-12,12), (13,-13), (-13,13)]

    fig, axs = plt.subplots(3, 3, sharey=True, sharex=True, figsize=(10, 7), layout='constrained')
    for dz, axss in zip(dzs, axs):
        for ihik, ax in zip(ihiks, axss):#itertools.chain.from_iterable(axs)):
            datadf = dfall.filter((pl.col("ih") == ihik[0]) & (pl.col("ik") == ihik[1]))
            print(datadf)
            for name in selmethods:
                if name == "orig+" or name == "orig":
                    methoddf = datadf.filter(pl.col("method")==name)
                else:
                    methoddf = datadf.filter((pl.col("method")==name) & (pl.col("dz")==dz))
                ax.plot(methoddf["angle"], methoddf["intensity"], linewidth=1, label=name, color=cmap[name][0], linestyle=cmap[name][2])
                if name != "orig+":
                    ax.scatter(methoddf["angle"], methoddf["intensity"], linewidth=1, s=20, label=name, color=cmap[name][0], marker=cmap[name][1])
            #ax.scatter(methoddf["angle"], methoddf["intensity"], color=cmap[name][0], marker=cmap[name][1])
            #methoddfop = datadf.filter(pl.col("method")=="orig+")
            #ax.plot(methoddf["angle"], abs(methoddf["intensity"]-methoddfop["intensity"]), label=name, color=cmap[name][0])
            #ax.scatter(methoddf["angle"], abs(methoddf["intensity"]-methoddfop["intensity"]), color=cmap[name][0], marker=cmap[name][1])
            ax.set_title("{} / {}".format(ihik,dz))
    #axs[0].text(18, 2e-4, "orig", color=cmap["orig"][0])
    #axs[0].text(1, 0.01, "opt", color=cmap["opt"][0])
    #axs[0].text(0.07, 0.001, "rk4", color=cmap["rk4"][0])
    #axs[0].text(0.3, 1e-7, "sp4", color=cmap["sp4"][0])
    #axs[0].text(0.03, 1e-7, "sp6", color=cmap["sp6"][0])
    #axs[1].text(100, 2e-4, "orig", color=cmap["orig"][0])
    #axs[1].text(5, 0.008, "opt", color=cmap["opt"][0])
    #axs[1].text(0.5, 0.001, "rk4", color=cmap["rk4"][0])
    #axs[1].text(1.5, 1e-7, "sp4", color=cmap["sp4"][0])
    #axs[1].text(0.15, 1e-7, "sp6", color=cmap["sp6"][0])
    #axs[2].text(1600, 2e-4, "orig", color=cmap["orig"][0])
    #axs[2].text(200, 0.008, "opt", color=cmap["opt"][0])
    #axs[2].text(20, 0.003, "rk4", color=cmap["rk4"][0])
    #axs[2].text(80, 5e-7, "sp4", color=cmap["sp4"][0])
    #axs[2].text(8, 5e-7, "sp6", color=cmap["sp6"][0])
    lhandles = [mlines.Line2D([], [], color=v[0], marker=v[1], label=k, linestyle=v[2]) for k, v in cmap.items()]
    #plt.yscale("log")
    plt.ylim(top=1.2e-2, bottom=-1e-3)
    #plt.ylim(bottom=1e-8)
    fig.supxlabel("glancing angle")
    fig.supylabel("intensity")
    plt.figlegend(handles=lhandles, loc="outside center right")
    #plt.figlegend(handles=lhandles, loc=(0.87, 0.2))
    #plt.tight_layout(pad=0.5, w_pad=0.1, h_pad=0.1, rect=[0., 0., 0.86, 1.])
    #plt.show()
    plt.savefig("rock.pdf", dpi=300)

    #t1 = dfmaxet.filter(pl.col("method") != "sp2")
    #to = t1.filter(pl.col("method") == "orig")
    #t2 = t1.join(to, on=["data"], how="left")\
    #        .filter(pl.col("err-acc") <= pl.col("err-acc_right"))\
    #        .filter(pl.col("sec.") == pl.min("sec.").over(["data", "method"]))\
    #        .with_columns((pl.col("sec._right")/pl.col("sec.")).alias("times"))\
    #        .sort(["method"]).select(["data", "method", "sec.", "dz", "err-acc", "times"])

    ##t2.sort(["data", "method"]).write_csv(file="scaledata.csv", float_precision=15)
    #file=open("table.tex", "w")
    #t2.sort(["data"]).to_pandas().to_latex(file, index=False,float_format="%.15e",columns=["data", "method", "dz", "err-acc", "sec.", "times"])
    #file.close()
    #print(t2.sort(["data", "method"]).write_csv(float_precision=15))
    #print(t3)

    #f = (pn.ggplot(t2, pn.aes("method", "times", color="method", fill="method"))
    #     + pn.geom_bar(stat="identity")
    #     + pn.facet_wrap(facets="data")
    #     #+ pn.scale_y_log10()
    #     + pn.scale_color_manual(cmap)
    #     + pn.scale_fill_manual(cmap)
    #     + pn.ylab("scaling")
    #     + pn.theme(
    #         legend_title=None))
    #f.draw()
    #plt.show()
    #f.save("plot-bar.pdf", width=12, height=4)

    #bar = px.bar(t2.to_pandas(), x="method", y="times", color="method", facet_col="data", log_y=True)
    #bar.show()
