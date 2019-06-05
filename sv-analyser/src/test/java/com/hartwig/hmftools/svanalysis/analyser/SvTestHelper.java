package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.DEFAULT_PROXIMITY_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvSampleAnalyser.setSvCopyNumberData;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SE_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SE_START;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.FusionDisruptionAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.LinkFinder;
import com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods;
import com.hartwig.hmftools.svanalysis.analysis.SvUtilities;
import com.hartwig.hmftools.svanalysis.cn.CnDataLoader;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.cn.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvaConfig;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

public class SvTestHelper
{
    public String SampleId;
    public List<SvVarData> AllVariants;
    public SvaConfig Config;
    public SvClusteringMethods ClusteringMethods;
    public ClusterAnalyser Analyser;
    public CnDataLoader CnDataLoader;
    public FusionDisruptionAnalyser FusionAnalyser;

    private int mNextVarId;

    public SvTestHelper()
    {
        Config = new SvaConfig(DEFAULT_PROXIMITY_DISTANCE);

        ClusteringMethods = new SvClusteringMethods(Config.ProximityDistance);
        Analyser = new ClusterAnalyser(Config, ClusteringMethods);
        CnDataLoader = new CnDataLoader("", "", null);
        Analyser.setCnDataLoader(CnDataLoader);

        Analyser.setRunValidationChecks(true);

        FusionAnalyser = null;

        SampleId = "TEST";
        AllVariants = Lists.newArrayList();

        Analyser.setSampleData(SampleId, AllVariants);
        mNextVarId = 0;

        Configurator.setRootLevel(Level.DEBUG);
    }

    public void initialiseFusions(SvGeneTranscriptCollection geneTranscriptCollection)
    {
        FusionAnalyser = new FusionDisruptionAnalyser();
        FusionAnalyser.initialise(null, "", Config, geneTranscriptCollection);
        FusionAnalyser.setHasValidConfigData(true);
    }

    public final String nextVarId() { return String.format("%d", mNextVarId++); }
    public void logVerbose(boolean toggle)
    {
        Config.LogVerbose = toggle;
        Analyser.getChainFinder().setLogVerbose(toggle);
        Analyser.getLinkFinder().setLogVerbose(toggle);
    }

    public void preClusteringInit()
    {
        ClusteringMethods.populateChromosomeBreakendMap(AllVariants);

        addCopyNumberData();

        ClusteringMethods.annotateNearestSvData();
        LinkFinder.findDeletionBridges(ClusteringMethods.getChrBreakendMap());
        ClusteringMethods.setSimpleVariantLengths(SampleId);
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
        ClusteringMethods.clearLOHBreakendData(SampleId);
        Analyser.getClusters().clear();
    }

    public void mergeOnProximity()
    {
        ClusteringMethods.clusterByProximity(Analyser.getClusters());
    }

    public final List<SvCluster> getClusters() { return Analyser.getClusters(); }

    public static SvVarData createSv(final String varId, final String chrStart, final String chrEnd,
            long posStart, long posEnd, int orientStart, int orientEnd, StructuralVariantType type, final String insertSeq)
    {
        return createTestSv(varId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, type,
                2, 2, 1, 1, 1, insertSeq);
    }

    // for convenience
    public static SvVarData createDel(final String varId, final String chromosome, long posStart, long posEnd)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, 1, -1, DEL,
                2, 2, 1, 1, 1, "");
    }

    public static SvVarData createIns(final String varId, final String chromosome, long posStart, long posEnd)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, 1, -1, INS,
                2, 2, 1, 1, 1, "");
    }

    public static SvVarData createDup(final String varId, final String chromosome, long posStart, long posEnd)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, -1, 1, DUP,
                3, 3, 1, 1, 1, "");
    }

    public static SvVarData createInv(final String varId, final String chromosome, long posStart, long posEnd, int orientation)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, orientation, orientation, INV,
                orientation == 1 ? 4 : 3, orientation == 1 ? 3 : 4, 1, 1, 1, "");
    }

    public static SvVarData createSgl(final String varId, final String chromosome, long position, int orientation, boolean isNoneSegment)
    {
        SvVarData var = createTestSv(varId, chromosome, "0", position, -1, orientation, -1, SGL,
                3, 0, 1, 0, 1, "");

        return var;
    }

    public static SvVarData createBnd(final String varId, final String chrStart, long posStart, int orientStart, final String chrEnd, long posEnd, int orientEnd)
    {
        SvVarData var = createTestSv(varId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, BND,
                3, 3, 1, 1, 1, "");

        return var;
    }

    public static SvVarData createTestSv(final String varId, final String chrStart, final String chrEnd,
            long posStart, long posEnd, int orientStart, int orientEnd, StructuralVariantType type, double ploidy)
    {
        // let the copy number test data routine take care of setting CN and CN change data
        return createTestSv(varId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, type,
                0, 0, ploidy, ploidy, ploidy, "");
    }

    public static SvVarData createTestSv(final String varId, final String chrStart, final String chrEnd,
            long posStart, long posEnd, int orientStart, int orientEnd, StructuralVariantType type,
            double cnStart, double cnEnd, double cnChgStart, double cnChgEnd, double ploidy, final String insertSeq)
    {
        StructuralVariantData svData =
                ImmutableStructuralVariantData.builder()
                        .id(Integer.parseInt(varId))
                        .startChromosome(chrStart)
                        .endChromosome(chrEnd)
                        .startPosition(posStart)
                        .endPosition(posEnd)
                        .startOrientation((byte)orientStart)
                        .endOrientation((byte)orientEnd)
                        .startHomologySequence("")
                        .endHomologySequence("")
                        .startAF(1.0)
                        .endAF(1.0)
                        .ploidy(ploidy)
                        .adjustedStartAF(1.0)
                        .adjustedEndAF(1.0)
                        .adjustedStartCopyNumber(cnStart)
                        .adjustedEndCopyNumber(cnEnd)
                        .adjustedStartCopyNumberChange(cnChgStart)
                        .adjustedEndCopyNumberChange(cnChgEnd)
                        .insertSequence(insertSeq)
                        .type(type)
                        .filter("PASS")
                        .imprecise(false)
                        .qualityScore(0.0)
                        .event("")
                        .startTumorVariantFragmentCount(0)
                        .startTumorReferenceFragmentCount(0)
                        .startNormalVariantFragmentCount(0)
                        .startNormalReferenceFragmentCount(0)
                        .endTumorVariantFragmentCount(0)
                        .endTumorReferenceFragmentCount(0)
                        .endNormalVariantFragmentCount(0)
                        .endNormalReferenceFragmentCount(0)
                        .startIntervalOffsetStart(0)
                        .startIntervalOffsetEnd(0)
                        .endIntervalOffsetStart(0)
                        .endIntervalOffsetEnd(0)
                        .inexactHomologyOffsetStart(0)
                        .inexactHomologyOffsetEnd(0)
                        .startLinkedBy("")
                        .endLinkedBy("")
                        .vcfId("")
                        .startRefContext("")
                        .endRefContext("")
                        .recovered(false)
                        .recoveryMethod("")
                        .recoveryFilter("")
                        .insertSequenceAlignments("")
                        .insertSequenceRepeatClass("")
                        .insertSequenceRepeatType("")
                        .insertSequenceRepeatOrientation((byte)0)
                        .insertSequenceRepeatCoverage(0.0)
                        .startAnchoringSupportDistance(0)
                        .endAnchoringSupportDistance(0)
                        .build();

        SvVarData var = new SvVarData(svData);

        String startArm = getChromosomalArm(var.chromosome(true), var.position(true));

        String endArm;
        if(!var.isNullBreakend())
            endArm = getChromosomalArm(var.chromosome(false), var.position(false));
        else
            endArm = CHROMOSOME_ARM_P;

        var.setChromosomalArms(startArm, endArm);

        // by default
        var.setPloidyRecalcData(var.getSvData().ploidy(), var.getSvData().ploidy());

        return var;
    }

    public void setDefaultPloidyCalcData()
    {
        AllVariants.stream().forEach(x -> x.setPloidyRecalcData(x.getSvData().ploidy(), x.getSvData().ploidy()));
    }

    private double calcActualBaf(double copyNumber, double nonDisruptedAP)
    {
        double disruptedPloidy  = copyNumber - nonDisruptedAP;
        return disruptedPloidy >= nonDisruptedAP ? disruptedPloidy / copyNumber : nonDisruptedAP / copyNumber;
    }

    public void addCopyNumberData()
    {
        // use SV breakend data to re-create the copy number segments
        // assume CN is 2 at the telomere and by default Actual BAF = 0.5

        // NOTE: positions adjusted for orientation are not done correctly
        Map<String, List<SvCNData>> chrCnDataMap = CnDataLoader.getChrCnDataMap();
        final Map<String, List<SvBreakend>> chrBreakendMap = ClusteringMethods.getChrBreakendMap();
        Map<Integer,SvCNData[]> svIdCnDataMap = CnDataLoader.getSvIdCnDataMap();

        chrCnDataMap.clear();
        svIdCnDataMap.clear();

        double nonDisruptedAP = 1;
        double currentCopyNumber = 2;
        double currentActualBaf = 0.5;

        int cnId = 0;
        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            List<SvBreakend> breakendList = entry.getValue();
            List<SvCNData> cnDataList = Lists.newArrayList();
            chrCnDataMap.put(chromosome, cnDataList);

            long centromerePosition = SvUtilities.getChromosomalArmLength(chromosome, CHROMOSOME_ARM_P);
            long chromosomeLength = SvUtilities.CHROMOSOME_LENGTHS.get(chromosome);

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvVarData var = breakend.getSV();
                double ploidy = var.ploidy();

                SvCNData cnData = null;

                if (i == 0)
                {
                    if(breakend.orientation() == 1)
                    {
                        currentCopyNumber = nonDisruptedAP + ploidy;
                    }

                    currentActualBaf = calcActualBaf(currentCopyNumber, nonDisruptedAP);

                    // add telomere segment at start, and centromere as soon as the breakend crosses the centromere
                    if(breakend.arm() == CHROMOSOME_ARM_Q)
                    {
                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, 0, centromerePosition,
                                currentCopyNumber, TELOMERE.toString(), CENTROMERE.toString(),
                                1, currentActualBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);

                        extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, breakend.position() - 1,
                                currentCopyNumber, CENTROMERE.toString(), var.type().toString(),
                                1, currentActualBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        if(breakend.orientation() == 1)
                        {
                            // copy number and actual BAF as expected
                        }
                        else
                        {
                            // no chromatid running to telomere on this end
                            currentCopyNumber = nonDisruptedAP;
                            currentActualBaf = 1;
                        }

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, 0, breakend.position() - 1,
                                currentCopyNumber, TELOMERE.toString(), var.type().toString(),
                                1, currentActualBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                }

                if(breakend.orientation() == 1)
                {
                    // copy number drop
                    currentCopyNumber -= ploidy;
                }
                else
                {
                    // no chromatid running to telomere on this end
                    currentCopyNumber += ploidy;
                }

                currentActualBaf = calcActualBaf(currentCopyNumber, nonDisruptedAP);

                if (i < breakendList.size() - 1)
                {
                    final SvBreakend nextBreakend = breakendList.get(i + 1);

                    if(breakend.arm() == CHROMOSOME_ARM_P && nextBreakend.arm() == CHROMOSOME_ARM_Q)
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), centromerePosition-1,
                                currentCopyNumber, var.type().toString(), CENTROMERE.toString(),
                                1, currentActualBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnDataList.add(cnData);

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, nextBreakend.position() - 1,
                                currentCopyNumber, CENTROMERE.toString(), nextBreakend.getSV().type().toString(),
                                1, currentActualBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), nextBreakend.position() - 1,
                                currentCopyNumber, var.type().toString(), nextBreakend.getSV().type().toString(),
                                1, currentActualBaf, 100);

                        cnData.setIndex(cnDataList.size());
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
                                1, currentActualBaf, 100);

                        cnData.setIndex(cnDataList.size());
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
                                1, currentActualBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnDataList.add(cnData);
                    }
                }

                SvCNData[] cnDataPair = svIdCnDataMap.get(var.id());

                if(cnDataPair == null)
                {
                    cnDataPair = new SvCNData[2];
                    svIdCnDataMap.put(var.dbId(), cnDataPair);
                }

                cnDataPair[breakend.usesStart() ? SE_START : SE_END] = cnData;

                // set copy number data back into the SV
                double beCopyNumber = breakend.orientation() == 1 ? currentCopyNumber + ploidy : currentCopyNumber;
                breakend.getSV().setCopyNumberData(breakend.usesStart(), beCopyNumber, ploidy);
            }
        }

        setSvCopyNumberData(
                AllVariants,
                CnDataLoader.getSampleSvPloidyCalcMap().get(SampleId),
                CnDataLoader.getSvIdCnDataMap(),
                CnDataLoader.getChrCnDataMap());

    }

}
