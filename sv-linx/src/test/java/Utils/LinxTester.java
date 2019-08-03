package Utils;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.populateChromosomeBreakendMap;
import static com.hartwig.hmftools.linx.analysis.SvSampleAnalyser.setSvCopyNumberData;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvaConstants.DEFAULT_PROXIMITY_DISTANCE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.linx.analysis.ClusterAnalyser;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.analysis.SvUtilities;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

public class LinxTester
{
    public String SampleId;
    public List<SvVarData> AllVariants;
    public LinxConfig Config;
    public ClusterAnalyser Analyser;
    public CnDataLoader CnDataLoader;
    public FusionDisruptionAnalyser FusionAnalyser;

    private int mNextVarId;

    // assume an A-allele which is unaffected by the SVs, and a B-allele which is
    private double mOtherAllelePloidy;
    private double mUndisruptedAllelePloidy; // the ploidy of the undisrupted B-allele

    private static final Logger LOGGER = LogManager.getLogger(LinxTester.class);

    public LinxTester()
    {
        Config = new LinxConfig(DEFAULT_PROXIMITY_DISTANCE);

        Analyser = new ClusterAnalyser(Config);
        CnDataLoader = new CnDataLoader( "", null);
        Analyser.setCnDataLoader(CnDataLoader);

        Analyser.setRunValidationChecks(true);

        FusionAnalyser = null;

        SampleId = "TEST";
        AllVariants = Lists.newArrayList();

        Analyser.setSampleData(SampleId, AllVariants);
        mNextVarId = 0;

        mOtherAllelePloidy = 1;
        mUndisruptedAllelePloidy = 0;

        Configurator.setRootLevel(Level.DEBUG);
    }

    public void initialiseFusions(SvGeneTranscriptCollection geneTranscriptCollection)
    {
        FusionAnalyser = new FusionDisruptionAnalyser();
        FusionAnalyser.initialise(null, "", Config, geneTranscriptCollection);
        FusionAnalyser.setHasValidConfigData(true);
    }

    public final int nextVarId() { return mNextVarId++; }
    public void logVerbose(boolean toggle)
    {
        Config.LogVerbose = toggle;
        Analyser.getChainFinder().setLogVerbose(toggle);
        Analyser.getLinkFinder().setLogVerbose(toggle);
    }

    public void preClusteringInit()
    {
        preClusteringInit(false);
    }

    public void preClusteringInit(boolean includePloidyCalcs)
    {
        // have to manually trigger breakend map creation since the CN data creation uses it
        Analyser.getState().reset();
        populateChromosomeBreakendMap(AllVariants, Analyser.getState());

        populateCopyNumberData(includePloidyCalcs);

        Analyser.preClusteringPreparation();
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

    public boolean hasClusterWithSVs(final List<SvVarData> svList)
    {
        return findClusterWithSVs(svList) != null;
    }

    public final SvCluster findClusterWithSVs(final List<SvVarData> svList)
    {
        for(final SvCluster cluster : Analyser.getClusters())
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

    public void setNonClusterAllelePloidies(double otherAllele, double undisruptedAllele)
    {
        mOtherAllelePloidy = otherAllele;
        mUndisruptedAllelePloidy = undisruptedAllele;
    }

    public void populateCopyNumberData(boolean includePloidyCalcs)
    {
        createCopyNumberData();

        if(includePloidyCalcs)
            CnDataLoader.calculateAdjustedPloidy(SampleId);

        setSvCopyNumberData(
                AllVariants,
                CnDataLoader.getSvPloidyCalcMap(),
                CnDataLoader.getSvIdCnDataMap(),
                CnDataLoader.getChrCnDataMap());
    }

    private double calcActualBaf(double copyNumber)
    {
        if(copyNumber == 0)
            return 0;

        double bAllelePloidy = max(copyNumber - mOtherAllelePloidy, 0);

        if(bAllelePloidy >= mOtherAllelePloidy)
            return bAllelePloidy / copyNumber;
        else
            return mOtherAllelePloidy / copyNumber;
    }

    private void createCopyNumberData()
    {
        // use SV breakend data to re-create the copy number segments

        Map<String, List<SvCNData>> chrCnDataMap = CnDataLoader.getChrCnDataMap();
        final Map<String, List<SvBreakend>> chrBreakendMap = Analyser.getState().getChrBreakendMap();
        Map<Integer,SvCNData[]> svIdCnDataMap = CnDataLoader.getSvIdCnDataMap();

        chrCnDataMap.clear();
        svIdCnDataMap.clear();

        int cnId = 0;
        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            List<SvBreakend> breakendList = entry.getValue();
            List<SvCNData> cnDataList = Lists.newArrayList();
            chrCnDataMap.put(chromosome, cnDataList);

            // work out the net copy number from all SVs going out to P-arm telomere for the correct starting copy number
            double netSvPloidy = max(breakendList.stream().mapToDouble(x -> x.ploidy() * x.orientation()).sum(), 0);

            double currentCopyNumber = mOtherAllelePloidy + mUndisruptedAllelePloidy + netSvPloidy;

            long centromerePosition = SvUtilities.getChromosomalArmLength(chromosome, CHROMOSOME_ARM_P);
            long chromosomeLength = SvUtilities.CHROMOSOME_LENGTHS.get(chromosome);

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final StructuralVariantData svData = breakend.getSV().getSvData();
                final SvVarData var = breakend.getSV();
                double ploidy = var.ploidy();

                double ploidyChange = -ploidy * breakend.orientation();

                SvCNData cnData = null;

                if (i == 0)
                {
                    // =IF(A15="DUP",g15,-G15)+MAX(G12,0)
                    if(breakend.getSV().type() == DUP && breakendList.get(i + 1).getSV() == breakend.getSV())
                    {
                        // starts with a DUP so don't treat the first breakend as a copy-number drop
                        currentCopyNumber += +ploidyChange;
                    }
                    else
                    {
                        currentCopyNumber += max(-ploidyChange, 0);
                    }

                    if(currentCopyNumber < 0)
                    {
                        LOGGER.error("invalid copy number({}) at telomere", currentCopyNumber);
                        return;
                    }

                    double actualBaf = calcActualBaf(currentCopyNumber);

                    // add telomere segment at start, and centromere as soon as the breakend crosses the centromere
                    if(breakend.arm() == CHROMOSOME_ARM_Q)
                    {
                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, 0, centromerePosition,
                                currentCopyNumber, TELOMERE.toString(), CENTROMERE.toString(),
                                1, actualBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);

                        extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, breakend.position() - 1,
                                currentCopyNumber, CENTROMERE.toString(), var.type().toString(),
                                1, actualBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, 0, breakend.position() - 1,
                                currentCopyNumber, TELOMERE.toString(), var.type().toString(),
                                1, actualBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                }

                // orientation determines copy number drop or gain
                currentCopyNumber += ploidyChange;

                if(currentCopyNumber < 0)
                {
                    LOGGER.error("invalid copy number({}) at breakend({})", currentCopyNumber, breakend);
                    return;
                }

                double actualBaf = calcActualBaf(currentCopyNumber);

                if (i < breakendList.size() - 1)
                {
                    final SvBreakend nextBreakend = breakendList.get(i + 1);
                    final StructuralVariantData nextSvData = nextBreakend.getSV().getSvData();

                    if(breakend.arm() == CHROMOSOME_ARM_P && nextBreakend.arm() == CHROMOSOME_ARM_Q)
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), centromerePosition-1,
                                currentCopyNumber, var.type().toString(), CENTROMERE.toString(),
                                1, actualBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, nextBreakend.position() - 1,
                                currentCopyNumber, CENTROMERE.toString(), nextBreakend.getSV().type().toString(),
                                1, actualBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), nextBreakend.position() - 1,
                                currentCopyNumber, var.type().toString(), nextBreakend.getSV().type().toString(),
                                1, actualBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);
                    }
                }
                else
                {
                    // last breakend runs out to the telomere
                    if(breakend.arm() == CHROMOSOME_ARM_P)
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), centromerePosition - 1,
                                currentCopyNumber,
                                var.type().toString(), CENTROMERE.toString(),
                                1, actualBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, chromosomeLength,
                                breakend.getCopyNumber(false),
                                CENTROMERE.toString(), TELOMERE.toString(),
                                1, 0.5, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), chromosomeLength,
                                currentCopyNumber, var.type().toString(), TELOMERE.toString(),
                                1, actualBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);
                    }
                }

                SvCNData[] cnDataPair = svIdCnDataMap.get(var.id());

                if(cnDataPair == null)
                {
                    cnDataPair = new SvCNData[2];
                    svIdCnDataMap.put(var.id(), cnDataPair);
                }

                cnDataPair[breakend.usesStart() ? SE_START : SE_END] = cnData;

                // set copy number data back into the SV
                double beCopyNumber = breakend.orientation() == 1 ? currentCopyNumber + ploidy : currentCopyNumber;
                breakend.getSV().setCopyNumberData(breakend.usesStart(), beCopyNumber, ploidy);
            }
        }
    }

}
