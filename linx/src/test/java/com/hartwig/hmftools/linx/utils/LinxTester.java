package com.hartwig.hmftools.linx.utils;

import static com.hartwig.hmftools.common.purple.Gender.MALE;
import static com.hartwig.hmftools.linx.analysis.AnnotationExtension.DOUBLE_MINUTES;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.linkSglMappedInferreds;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.populateChromosomeBreakendMap;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.initialiseSV;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.analysis.AnnotationExtension;
import com.hartwig.hmftools.linx.analysis.ClusterAnalyser;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.annotators.LineElementAnnotator;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.CnSegmentBuilder;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.drivers.DriverGeneAnnotator;
import com.hartwig.hmftools.linx.fusion.FusionConfig;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.fusion.FusionResources;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisSampleData;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

public class LinxTester
{
    public String SampleId;
    public List<SvVarData> AllVariants;
    public LinxConfig Config;
    public ClusterAnalyser Analyser;
    public CnDataLoader CnDataLoader;
    public LineElementAnnotator LineAnnotator;
    public FusionDisruptionAnalyser FusionAnalyser;
    public DriverGeneAnnotator DriverAnnotator;
    public CohortDataWriter CohortWriter;
    public VisSampleData VisData;

    private CnSegmentBuilder mCnSegmentBuilder;
    private int mNextVarId;

    public LinxTester()
    {
        this(false);
    }

    public LinxTester(boolean isGermline)
    {
        Config = new LinxConfig(isGermline);
        LinxConfig.addConfig(Config.CmdLineConfig);
        FusionConfig.addConfig(Config.CmdLineConfig);
        Config.AnnotationExtensions.add(DOUBLE_MINUTES);

        Analyser = new ClusterAnalyser(Config, null);

        if(!isGermline)
        {
            CnDataLoader = new CnDataLoader("");
            Analyser.setCnDataLoader(CnDataLoader);

            mCnSegmentBuilder = new CnSegmentBuilder();
        }
        else
        {
            mCnSegmentBuilder = null;
        }

        SampleId = "TEST";
        VisData = new VisSampleData();

        CohortWriter = new CohortDataWriter(Config, null);

        LineAnnotator = new LineElementAnnotator(Config.ProximityDistance);
        Analyser.setLineAnnotator(LineAnnotator);

        Analyser.setRunValidationChecks(true);

        FusionAnalyser = null;

        AllVariants = Lists.newArrayList();

        Analyser.setSampleData(SampleId, AllVariants);
        mNextVarId = 0;

        // Configurator.setRootLevel(Level.DEBUG);
    }

    public void initialiseFusions(final EnsemblDataCache ensemblDataCache)
    {
        FusionAnalyser = new FusionDisruptionAnalyser(
                Config, ensemblDataCache, new FusionResources(Config.CmdLineConfig), CohortWriter, VisData);
    }

    public void initialiseDriverGeneAnnotator(final EnsemblDataCache ensemblDataCache)
    {
        DriverAnnotator = new DriverGeneAnnotator(ensemblDataCache, Config, CnDataLoader, CohortWriter, VisData);
    }

    public final int nextVarId() { return mNextVarId++; }
    public void logVerbose(boolean toggle)
    {
        Config.LogVerbose = toggle;

        if(toggle)
            Configurator.setRootLevel(Level.TRACE);

        Analyser.getChainFinder().setLogVerbose(toggle);
    }

    public void addAndCluster(SvVarData var1, SvVarData var2)
    {
        clearClustersAndSVs();
        AllVariants.add(var1);
        AllVariants.add(var2);
        preClusteringInit();
        Analyser.clusterAndAnalyse();
    }

    public void preClusteringInit()
    {
        preClusteringInit(false);
    }

    public void preClusteringInit(boolean includePloidyCalcs)
    {
        // have to manually trigger breakend map creation since the CN data creation uses it
        Analyser.getState().reset();
        linkSglMappedInferreds(AllVariants);
        AllVariants.stream().filter(x -> x.getLinkedSVs() != null).filter(x -> !x.isSglBreakend()).forEach(x -> initialiseSV(x));

        populateChromosomeBreakendMap(AllVariants, Analyser.getState());

        if(!Config.IsGermline)
            populateCopyNumberData(includePloidyCalcs);

        Analyser.setSampleData(SampleId, AllVariants);
        Analyser.preClusteringPreparation();
    }

    public void addLohEvent(final SvBreakend breakend1, final SvBreakend breakend2)
    {
        if(breakend1.orientation() != 1 || breakend2.orientation() != -1 || !breakend1.chromosome().equals(breakend2.chromosome()))
            return;

        LohEvent lohEvent = new LohEvent(breakend1.chromosome(), breakend1.position(), breakend2.position(),
                breakend1.getSV().typeStr(), breakend2.getSV().typeStr(), 1, breakend1.getSV().id(), breakend2.getSV().id());

        CnDataLoader.getLohData().add(lohEvent);
    }

    public void addClusterAndSVs(final SvCluster cluster)
    {
        Analyser.getClusters().add(cluster);
        AllVariants.addAll(cluster.getSVs());
    }

    public void clearClustersAndSVs()
    {
        // in case SVs are to be used again and re-clustered
        for(SvVarData var : AllVariants)
        {
            var.setCluster(null);
        }

        AllVariants.clear();
        Analyser.getState().reset();
        Analyser.getClusters().clear();
    }

    public final List<SvCluster> getClusters() { return Analyser.getClusters(); }
    public final List<SvCluster> getAllClusters() { return Analyser.getAllClusters(); }

    public boolean hasClusterWithSVs(final List<SvVarData> svList)
    {
        return findClusterWithSVs(svList) != null;
    }

    public final SvCluster findClusterWithSVs(final List<SvVarData> svList)
    {
        for(final SvCluster cluster : Analyser.getAllClusters())
        {
            if(cluster.getSvCount() != svList.size())
                continue;

            boolean hasAll = true;

            for(final SvVarData var : svList)
            {
                if(!cluster.getSVs().contains(var))
                {
                    hasAll = false;
                    break;
                }
            }

            if(hasAll)
                return cluster;
        }

        return null;
    }

    public final SvChain findChainWithSVs(final SvCluster cluster, final List<SvVarData> svList)
    {
        for(final SvChain chain : cluster.getChains())
        {
            if(chain.getSvList().size() != svList.size())
                continue;

            boolean hasAll = true;

            for(final SvVarData var : svList)
            {
                if(!chain.getSvList().contains(var))
                {
                    hasAll = false;
                    break;
                }
            }

            if(hasAll)
                return chain;
        }

        return null;
    }

    public void setNonClusterAllelePloidies(double otherAllele, double undisruptedAllele)
    {
        mCnSegmentBuilder.setAllelePloidies(otherAllele, undisruptedAllele);
    }

    public void populateCopyNumberData(boolean includePloidyCalcs)
    {
        mCnSegmentBuilder.createCopyNumberData(CnDataLoader, Analyser.getState().getChrBreakendMap());

        if(includePloidyCalcs)
            CnDataLoader.calculateAdjustedJcn(SampleId);

        mCnSegmentBuilder.setSamplePurity(CnDataLoader, 1, 2, MALE);

        CnDataLoader.createChrCopyNumberMap();

        SvCNData.setSvCopyNumberData(
                AllVariants,
                CnDataLoader.getSvJcnCalcMap(),
                CnDataLoader.getSvIdCnDataMap(),
                CnDataLoader.getChrCnDataMap());
    }

}
