package com.hartwig.hmftools.amber;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.LOH_CENTROMERE_LENGTH;
import static com.hartwig.hmftools.amber.AmberConstants.LOH_MIN_LENGTH;
import static com.hartwig.hmftools.amber.AmberConstants.LOH_MOVING_AVG_COUNT;
import static com.hartwig.hmftools.amber.AmberConstants.LOH_MOVING_ITEM_COUNT;
import static com.hartwig.hmftools.amber.AmberConstants.LOH_POINT_DENSITY;
import static com.hartwig.hmftools.amber.AmberUtils.aboveQualFilter;
import static com.hartwig.hmftools.common.amber.AmberGender.determineGender;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.Gender;

public class LohIdentifier
{
    private final RefGenomeCoordinates mGenomeCoordinates;
    private final int mTumorMinDepth;

    public LohIdentifier(final RefGenomeVersion refGenVersion, final int tumorMinDepth)
    {
        mGenomeCoordinates = RefGenomeCoordinates.refGenomeCoordinates(refGenVersion);
        mTumorMinDepth = tumorMinDepth;
    }

    public void findLohSegments(final ListMultimap<Chromosome, TumorBAF> chromosomeBafs)
    {
        for(Chromosome chromosome : chromosomeBafs.keySet())
        {
            // DISCUSS
            //if(chromosome == HumanChromosome._Y)
            //    continue;

            if(chromosome.isAllosome())
                continue;

            List<TumorBAF> tumorBAFS = chromosomeBafs.get(chromosome);
            findChromosomeLohSegments(chromosome, tumorBAFS);
        }
    }

    private void findChromosomeLohSegments(final Chromosome chromosome, final List<TumorBAF> bafs)
    {
        List<PositionEvidence> positionEvidences = bafs.stream()
                .filter(x -> x.TumorEvidence.ReadDepth >= mTumorMinDepth)
                .filter(x -> aboveQualFilter(x.TumorEvidence))
                .map(x -> x.TumorEvidence).collect(Collectors.toList());

        for(PositionEvidence posEvidence : positionEvidences)
        {
            posEvidence.setLohData();
        }

        // calculate a moving average across all points
        int validCountHigh = 0;
        int validCountMed = 0;
        int validCountLow = 0;

        for(int i = 0; i < positionEvidences.size(); ++i)
        {
            PositionEvidence posEvidence = positionEvidences.get(i);

            if(posEvidence.LohData.ValidBafHigh)
                ++validCountHigh;

            if(posEvidence.LohData.ValidBafMedium)
                ++validCountMed;

            if(posEvidence.LohData.ValidBafLow)
                ++validCountLow;

            if(i >= LOH_MOVING_ITEM_COUNT)
            {
                PositionEvidence purgedTumorBAF = positionEvidences.get(i - LOH_MOVING_ITEM_COUNT);

                if(purgedTumorBAF.LohData.ValidBafHigh)
                    --validCountHigh;

                if(purgedTumorBAF.LohData.ValidBafMedium)
                    --validCountMed;

                if(purgedTumorBAF.LohData.ValidBafLow)
                    --validCountLow;

                if(validCountHigh < 0 || validCountMed < 0 || validCountLow < 0)
                {
                    AMB_LOGGER.error("chromosome({}) invalid valid BAF counts", chromosome);
                }

                PositionEvidence targetTumorBAF = positionEvidences.get(i - LOH_MOVING_AVG_COUNT);
                targetTumorBAF.LohData.setValidBafData(validCountHigh, validCountMed, validCountLow);
            }
        }

        // find continuous regions of LOH
        int centromere = mGenomeCoordinates.Centromeres.get(chromosome);

        int lohIndexStart = -1;
        int lohIndexEnd = -1;

        for(int i = 0; i < positionEvidences.size(); ++i)
        {
            PositionEvidence posEvidence = positionEvidences.get(i);

            if(posEvidence.LohData.LohCandidate)
            {
                if(lohIndexStart >= 0)
                {
                    lohIndexEnd = i;
                }
                else
                {
                    lohIndexStart = lohIndexEnd = i;
                }
            }
            else if(lohIndexStart >= 0)
            {
                PositionEvidence lohBafStart = positionEvidences.get(lohIndexStart);
                PositionEvidence lohBafEnd = positionEvidences.get(lohIndexEnd);

                int lohLength = lohBafEnd.position() - lohBafStart.position();

                if(centromere > lohBafStart.position() && centromere < lohBafEnd.position())
                {
                    lohLength -= LOH_CENTROMERE_LENGTH;
                }

                if(lohLength >= LOH_MIN_LENGTH)
                {
                    AMB_LOGGER.debug("chromosome({}) LOH segment({} - {}) points({})",
                            chromosome, lohBafStart.position(), lohBafEnd.position(), lohIndexEnd - lohIndexStart + 1);

                    // mark these points as LOH, plus the ones either side used in the moving average
                    int widenedIndexStart = max(lohIndexStart - LOH_MOVING_AVG_COUNT, 0);
                    int widenedIndexEnd = min(lohIndexEnd + LOH_MOVING_AVG_COUNT, positionEvidences.size() - 1);

                    int lastRecoveredPoint = -1;

                    for(int j = widenedIndexStart; j <= widenedIndexEnd; ++j)
                    {
                        PositionEvidence lohPosEvidence = positionEvidences.get(j);

                        if(lohPosEvidence.ReadDepth >= mTumorMinDepth)
                        {
                            if(lastRecoveredPoint < 0 || j - lastRecoveredPoint >= LOH_POINT_DENSITY)
                            {
                                lohPosEvidence.LohData.LohRecovered = true;
                                lastRecoveredPoint = j;
                            }
                        }
                    }
                }

                lohIndexStart = lohIndexEnd = -1; // reset
            }
        }
    }

    public static List<AmberBAF> checkGenderSites(final RefGenomeVersion refGenomeVersion, final List<AmberBAF> amberBAFList)
    {
        ListMultimap<Chromosome, AmberBAF> chrBafMap = ArrayListMultimap.create();
        amberBAFList.forEach(x -> chrBafMap.put(HumanChromosome.fromString(x.Chromosome), x));

        Gender gender = determineGender(refGenomeVersion, chrBafMap, null);

        if(gender == Gender.FEMALE)
            return amberBAFList;

        // remove the recovered sites on X
        return amberBAFList.stream().filter(x -> !(x.Chromosome.contains("X") && x.LohRecovered)).collect(Collectors.toList());
    }
}
