import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

with pl.StringCache():
    methodnames = {"99":"orig", "0":"opt", "7":"rk4", "4":"sp2", "5":"sp4", "6":"sp6"}
    datanames = ["test23", "test47", "test521"]#, "test1229"]
    datanames2 = ["n="+t[4:] for t in datanames]
    pl.Series(methodnames.values()).cast(pl.Categorical)
    pl.Series(datanames2).cast(pl.Categorical)
    dfer = pl.DataFrame()
    for t in datanames:
        d = pl.read_csv(
                file = "../"+t+"/test-results/perf-comp-orig.txt",
                has_header = False,
                columns = [2,3,4,5,6],
                sep = " ",
                dtypes = [pl.Categorical, pl.Float64, pl.Float64, pl.Float64, pl.Float64])
        d.columns = ["method", "dz", "angle", "err-orig", "erel-orig"]
        d = d.with_columns(pl.lit("n="+t[4:]).alias("data")).with_columns(pl.col("method").apply(lambda x: methodnames[x]).cast(pl.Categorical))
        dfer = dfer.vstack(d)

    dfac = pl.DataFrame()
    for t in datanames:
        d = pl.read_csv(
                file = "../"+t+"/test-results/perf-comp-acc.txt",
                has_header = False,
                columns = [2,3,4,5,6],
                sep = " ",
                dtypes = [pl.Categorical, pl.Float64, pl.Float64, pl.Float64, pl.Float64])
        d.columns = ["method", "dz", "angle", "err-acc", "erel-acc"]
        d = d.with_columns(pl.lit("n="+t[4:]).alias("data")).with_columns(pl.col("method").apply(lambda x: methodnames[x]).cast(pl.Categorical))
        dfac = dfac.vstack(d)
    dfer = dfer.join(dfac, on=["method", "dz", "angle", "data"], how="outer")

    dftm = pl.DataFrame()
    for t in datanames:
        d = pl.read_csv(
                file = "../"+t+"/test-results/time.txt",
                has_header = False,
                sep = " ",
                dtypes = [pl.Categorical, pl.Float64, pl.Int32, pl.Float64])
        d.columns = ["method", "dz", "trial", "sec."]
        d = d.with_columns(pl.lit("n="+t[4:]).alias("data")).with_columns(pl.col("method").apply(lambda x: methodnames[x]).cast(pl.Categorical))
        dftm = dftm.vstack(d)
    dftmag = dftm.groupby(["method", "dz", "data"]).agg([pl.mean("sec.").alias("sec."), pl.std("sec.").alias("sec. std")])

    dfmaxe = dfer.groupby(["method", "dz", "data"], maintain_order=True).max()
    print(dfmaxe)

    dfmaxet = dfmaxe.join(dftmag, on=["method", "dz", "data"]).sort(["data", "method", "dz"])
    #t0 = dfmaxet.filter(pl.col("method") != "sp2")
    t0 = dftm.join(dfmaxe, on=["method", "dz", "data"]).filter(pl.col("method") != "sp2")
    t0mean = t0.groupby(["method", "dz", "data"], maintain_order=True).median()

    pltcmap = lambda x: "#{:02x}{:02x}{:02x}".format(*tuple(((np.array(plt.get_cmap("Set1").colors)*255)[x].astype(int))))
    cmap = {
            "orig":(pltcmap(0),"o"),
            "opt":(pltcmap(1),"."),
            "rk4":(pltcmap(2),"*"),
            "sp4":(pltcmap(3),"+"),
            "sp6":(pltcmap(4),"x")}

    fig, axs = plt.subplots(1, 3, sharey=True, figsize=(8, 3))
    for data, ax in zip(datanames2, axs):
        datadf = t0mean.filter(pl.col("data") == data)
        for name, methoddf in datadf.groupby("method"):
            sc = t0.filter((pl.col("data")==data) & (pl.col("method") == name))
            ax.plot(methoddf["sec."], methoddf["err-acc"], label=name, color=cmap[name][0])
            ax.scatter(sc["sec."], sc["err-acc"], color=cmap[name][0], marker=cmap[name][1])
        ax.set_xscale("log")
        ax.set_title(data)
    axs[0].text(20, 0.002, "orig", color=cmap["orig"][0])
    axs[0].text(1, 0.01, "opt", color=cmap["opt"][0])
    axs[0].text(0.15, 0.0001, "rk4", color=cmap["rk4"][0])
    axs[0].text(0.3, 0.5e-8, "sp4", color=cmap["sp4"][0])
    axs[0].text(0.03, 0.5e-9, "sp6", color=cmap["sp6"][0])
    axs[1].text(120, 0.002, "orig", color=cmap["orig"][0])
    axs[1].text(5, 0.008, "opt", color=cmap["opt"][0])
    axs[1].text(0.9, 0.0001, "rk4", color=cmap["rk4"][0])
    axs[1].text(1.5, 0.5e-8, "sp4", color=cmap["sp4"][0])
    axs[1].text(0.15, 0.1e-9, "sp6", color=cmap["sp6"][0])
    axs[2].text(1800, 0.002, "orig", color=cmap["orig"][0])
    axs[2].text(200, 0.008, "opt", color=cmap["opt"][0])
    axs[2].text(40, 0.0003, "rk4", color=cmap["rk4"][0])
    axs[2].text(80, 0.5e-8, "sp4", color=cmap["sp4"][0])
    axs[2].text(8, 0.1e-9, "sp6", color=cmap["sp6"][0])
    lhandles = [mlines.Line2D([], [], color=v[0], marker=v[1], label=k) for k, v in cmap.items()]
    plt.yscale("log")
    fig.supxlabel("time in sec.")
    fig.supylabel("eacc")
    plt.figlegend(handles=lhandles, loc="center right")
    plt.tight_layout(pad=0.5, w_pad=0.1, h_pad=0.1, rect=[0., 0., 0.91, 1.])
    #plt.show()
    plt.savefig("plot-acc.pdf", dpi=300)


    fig, axs = plt.subplots(1, 3, sharey=True, figsize=(8, 3))
    for data, ax in zip(datanames2, axs):
        datadf = t0mean.filter(pl.col("data") == data)
        for name, methoddf in datadf.groupby("method"):
            sc = t0.filter((pl.col("data")==data) & (pl.col("method") == name))
            if name=="orig":
                # XXX eorig of orig is exactly zero. replace value to small one (and write so in article)
                ax.plot(methoddf["sec."], [1e-4]*(methoddf["sec."].len()), label=name, color=cmap[name][0])
                ax.scatter(sc["sec."], [1e-4]*len(sc["sec."]), color=cmap[name][0], marker=cmap[name][1])
            else:
                ax.plot(methoddf["sec."], methoddf["err-orig"], label=name, color=cmap[name][0])
                ax.scatter(sc["sec."], sc["err-orig"], color=cmap[name][0], marker=cmap[name][1])
        ax.set_xscale("log")
        ax.set_title(data)
    axs[0].text(18, 1.5e-4, "orig", color=cmap["orig"][0])
    axs[0].text(1, 0.01, "opt", color=cmap["opt"][0])
    axs[0].text(0.07, 0.001, "rk4", color=cmap["rk4"][0])
    axs[0].text(0.3, 2e-4, "sp4", color=cmap["sp4"][0])
    axs[0].text(0.03, 2e-4, "sp6", color=cmap["sp6"][0])
    axs[1].text(100, 1.5e-4, "orig", color=cmap["orig"][0])
    axs[1].text(5, 0.008, "opt", color=cmap["opt"][0])
    axs[1].text(0.5, 0.001, "rk4", color=cmap["rk4"][0])
    axs[1].text(1.5, 1.5e-4, "sp4", color=cmap["sp4"][0])
    axs[1].text(0.15, 1.5e-4, "sp6", color=cmap["sp6"][0])
    axs[2].text(1600, 1.5e-4, "orig", color=cmap["orig"][0])
    axs[2].text(200, 0.008, "opt", color=cmap["opt"][0])
    axs[2].text(20, 0.003, "rk4", color=cmap["rk4"][0])
    axs[2].text(80, 2.5e-4, "sp4", color=cmap["sp4"][0])
    axs[2].text(8, 2.5e-4, "sp6", color=cmap["sp6"][0])
    lhandles = [mlines.Line2D([], [], color=v[0], marker=v[1], label=k) for k, v in cmap.items()]
    plt.yscale("log")
    plt.ylim(bottom=1e-4)
    fig.supxlabel("time in sec.")
    fig.supylabel("eorig")
    plt.figlegend(handles=lhandles, loc="center right")
    plt.tight_layout(pad=0.5, w_pad=0.1, h_pad=0.1, rect=[0., 0., 0.91, 1.])
    #plt.show()
    plt.savefig("plot-orig.pdf", dpi=300)

    t1 = dfmaxet.filter(pl.col("method") != "sp2")
    to = t1.filter(pl.col("method") == "orig")
    t2 = t1.join(to, on=["data"], how="left")\
            .filter(pl.col("err-acc") <= pl.col("err-acc_right"))\
            .filter(pl.col("sec.") == pl.min("sec.").over(["data", "method"]))\
            .with_columns((pl.col("sec._right")/pl.col("sec.")).alias("times"))\
            .sort(["method"]).select(["data", "method", "sec.", "dz", "err-acc", "times"])

    #t2.sort(["data", "method"]).write_csv(file="scaledata.csv", float_precision=15)
    file=open("table.tex", "w")
    t2.sort(["data"]).to_pandas().to_latex(file, index=False,float_format="%.15e",columns=["data", "method", "dz", "err-acc", "sec.", "times"])
    file.close()
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
