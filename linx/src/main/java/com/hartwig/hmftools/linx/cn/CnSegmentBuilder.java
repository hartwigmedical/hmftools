package com.hartwig.hmftools.linx.cn;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.purple.purity.FittedPurityMethod.NORMAL;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus.UNKNOWN;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.Q_ARM;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.linx.analysis.SvUtilities;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class CnSegmentBuilder
{
    // assume an A-allele which is unaffected by the SVs, and a B-allele which is
    private double mOtherAlleleJcn;
    private double mUndisruptedAlleleJcn; // the JCN of the undisrupted B-allele

    public CnSegmentBuilder()
    {
        mOtherAlleleJcn = 1;
        mUndisruptedAlleleJcn = 0;
    }

    public void setAllelePloidies(double otherAllele, double undisruptedAllele)
    {
        mOtherAlleleJcn = otherAllele;
        mUndisruptedAlleleJcn = undisruptedAllele;
    }

    public void createCopyNumberData(final CnDataLoader cnDataLoader, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        // use SV breakend data to re-create the copy number segments

        final Map<String, List<SvCNData>> chrCnDataMap = cnDataLoader.getChrCnDataMap();
        final Map<Integer,SvCNData[]> svIdCnDataMap = cnDataLoader.getSvIdCnDataMap();

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
            double netSvJcn = max(breakendList.stream().mapToDouble(x -> x.jcn() * x.orientation()).sum(), 0);

            double currentCopyNumber = mOtherAlleleJcn + mUndisruptedAlleleJcn + netSvJcn;

            int centromerePosition = SvUtilities.getChromosomalArmLength(chromosome, P_ARM);
            int chromosomeLength = SvUtilities.getChromosomeLength(chromosome);

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final StructuralVariantData svData = breakend.getSV().getSvData();
                final SvVarData var = breakend.getSV();
                double jcn = var.jcn();

                double jcnChange = -jcn * breakend.orientation();

                SvCNData cnData = null;

                if (i == 0)
                {
                    if(breakend.type() == DUP && breakendList.get(i + 1).getSV() == breakend.getSV())
                    {
                        // starts with a DUP so don't treat the first breakend as a copy-number drop
                        currentCopyNumber += jcnChange;
                    }
                    else
                    {
                        currentCopyNumber += max(-jcnChange, 0);
                    }

                    if(currentCopyNumber < 0)
                    {
                        LNX_LOGGER.error("invalid copy number({}) at telomere", currentCopyNumber);
                        return;
                    }

                    double actualBaf = calcActualBaf(currentCopyNumber);

                    // add telomere segment at start, and centromere as soon as the breakend crosses the centromere
                    if(breakend.arm() == Q_ARM)
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
                currentCopyNumber += jcnChange;

                if(currentCopyNumber < 0)
                {
                    LNX_LOGGER.warn("invalid copy number({}) from jcnChange({}) at breakend({})",
                            currentCopyNumber, formatJcn(jcnChange), breakend);
                    currentCopyNumber = 0;
                }

                double actualBaf = calcActualBaf(currentCopyNumber);

                if (i < breakendList.size() - 1)
                {
                    final SvBreakend nextBreakend = breakendList.get(i + 1);

                    if(breakend.arm() == P_ARM && nextBreakend.arm() == Q_ARM)
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), centromerePosition-1,
                                currentCopyNumber, var.type().toString(), CENTROMERE.toString(),
                                1, actualBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, nextBreakend.position() - 1,
                                currentCopyNumber, CENTROMERE.toString(), nextBreakend.type().toString(),
                                1, actualBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), nextBreakend.position() - 1,
                                currentCopyNumber, var.type().toString(), nextBreakend.type().toString(),
                                1, actualBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);
                    }
                }
                else
                {
                    // last breakend runs out to the telomere
                    if(breakend.arm() == P_ARM)
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), centromerePosition - 1,
                                currentCopyNumber,
                                var.type().toString(), CENTROMERE.toString(),
                                1, actualBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, chromosomeLength,
                                currentCopyNumber,
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
                double beCopyNumber = breakend.orientation() == 1 ? currentCopyNumber + jcn : currentCopyNumber;
                breakend.getSV().setCopyNumberData(breakend.usesStart(), beCopyNumber, jcn);
            }
        }
    }

    private double calcActualBaf(double copyNumber)
    {
        if(copyNumber == 0)
            return 0;

        double bAlleleJcn = max(copyNumber - mOtherAlleleJcn, 0);

        if(bAlleleJcn >= mOtherAlleleJcn)
            return bAlleleJcn / copyNumber;
        else
            return mOtherAlleleJcn / copyNumber;
    }

    public void createIndependentCopyNumberData(final CnDataLoader cnDataLoader, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        // set copy number data for each breakend irrespective of the breakends around it

        final Map<String, List<SvCNData>> chrCnDataMap = cnDataLoader.getChrCnDataMap();
        final Map<Integer,SvCNData[]> svIdCnDataMap = cnDataLoader.getSvIdCnDataMap();

        chrCnDataMap.clear();
        svIdCnDataMap.clear();

        double assumedCn = 2;
        double assumedBaf = 0.5;

        int cnId = 0;
        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();
            final List<SvCNData> cnDataList = Lists.newArrayList();
            chrCnDataMap.put(chromosome, cnDataList);

            int centromerePosition = SvUtilities.getChromosomalArmLength(chromosome, P_ARM);
            int chromosomeLength = SvUtilities.getChromosomeLength(chromosome);

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final StructuralVariantData svData = breakend.getSV().getSvData();
                final SvVarData var = breakend.getSV();
                double jcn = var.jcn();

                double jcnChange = -jcn * breakend.orientation();

                SvCNData cnData = null;

                if (i == 0)
                {
                    // add telomere segment at start, and centromere as soon as the breakend crosses the centromere
                    if(breakend.arm() == Q_ARM)
                    {
                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, 0, centromerePosition,
                                assumedCn, TELOMERE.toString(), CENTROMERE.toString(), 1, assumedBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);

                        extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, breakend.position() - 1,
                                assumedCn, CENTROMERE.toString(), var.type().toString(), 1, assumedBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, 0, breakend.position() - 1,
                                assumedCn, TELOMERE.toString(), var.type().toString(), 1, assumedBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                }

                if (i < breakendList.size() - 1)
                {
                    final SvBreakend nextBreakend = breakendList.get(i + 1);

                    if(breakend.arm() == P_ARM && nextBreakend.arm() == Q_ARM)
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), centromerePosition-1,
                                assumedCn, var.type().toString(), CENTROMERE.toString(), 1, assumedBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, nextBreakend.position() - 1,
                                assumedCn, CENTROMERE.toString(), nextBreakend.type().toString(), 1, assumedBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), nextBreakend.position() - 1,
                                assumedCn, var.type().toString(), nextBreakend.type().toString(), 1, assumedBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);
                    }
                }
                else
                {
                    // last breakend runs out to the telomere
                    if(breakend.arm() == P_ARM)
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), centromerePosition - 1,
                                assumedCn, var.type().toString(), CENTROMERE.toString(), 1, assumedBaf, 100);

                        cnData.setIndex(cnDataList.size());
                        cnData.setStructuralVariantData(svData, breakend.usesStart());
                        cnDataList.add(cnData);

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, chromosomeLength,
                                assumedCn, CENTROMERE.toString(), TELOMERE.toString(), 1, assumedBaf, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), chromosomeLength,
                                assumedCn, var.type().toString(), TELOMERE.toString(), 1, assumedBaf, 100);

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
                breakend.getSV().setCopyNumberData(breakend.usesStart(), assumedCn, jcn);
            }
        }
    }


    public void setSamplePurity(final CnDataLoader cnDataLoader, double purity, double ploidy, Gender gender)
    {
        FittedPurity fittedPurity = ImmutableFittedPurity.builder()
                .purity(purity)
                .ploidy(ploidy)
                .diploidProportion(1)
                .normFactor(1)
                .score(1)
                .somaticPenalty(0)
                .build();

        FittedPurityScore purityScore = ImmutableFittedPurityScore.builder()
                .maxPurity(1)
                .minPurity(1)
                .maxDiploidProportion(1)
                .maxPloidy(2)
                .minDiploidProportion(0)
                .minPloidy(2)
                .build();

        PurpleQC qc = ImmutablePurpleQC.builder()
                .amberGender(gender)
                .cobaltGender(gender)
                .purity(purity)
                .contamination(0)
                .deletedGenes(0)
                .copyNumberSegments(0)
                .unsupportedCopyNumberSegments(0)
                .method(NORMAL)
                .amberMeanDepth(0)
                .build();

        PurityContext purityContext = ImmutablePurityContext.builder()
                .bestFit(fittedPurity)
                .gender(gender)
                .microsatelliteIndelsPerMb(0)
                .microsatelliteStatus(UNKNOWN)
                .polyClonalProportion(0)
                .score(purityScore)
                .method(NORMAL)
                .version("1.0")
                .wholeGenomeDuplication(false)
                .tumorMutationalLoad(0)
                .tumorMutationalLoadStatus(TumorMutationalStatus.UNKNOWN)
                .tumorMutationalBurdenPerMb(0)
                .tumorMutationalBurdenStatus(TumorMutationalStatus.UNKNOWN)
                .svTumorMutationalBurden(0)
                .qc(qc)
                .build();

        cnDataLoader.setPurityContext(purityContext);
    }


}
