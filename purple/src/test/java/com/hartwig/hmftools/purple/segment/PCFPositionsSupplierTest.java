package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.utils.pcf.PCFSource.REFERENCE_RATIO;
import static com.hartwig.hmftools.common.utils.pcf.PCFSource.TUMOR_RATIO;
import static com.hartwig.hmftools.purple.FittingTestUtils.buildCobaltChromosomes;
import static com.hartwig.hmftools.purple.FittingTestUtils.createCobaltRatio;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;

import org.junit.Before;
import org.junit.Test;

public class PCFPositionsSupplierTest
{
    private final Chromosome chr = HumanChromosome._21;
    private AmberData mAmberData;
    private CobaltData mCobaltData;

    @Before
    public void init()
    {
        mAmberData = new AmberData(100, Gender.MALE);

        CobaltChromosomes cobaltChromosomes = buildCobaltChromosomes();
        mCobaltData = new CobaltData(cobaltChromosomes);
    }

    @Test
    public void amberAndCobaltBothEmpty()
    {
        assertTrue(PCFPositionsSupplier.createPositions(mAmberData, mCobaltData).isEmpty());
    }

    @Test
    public void amberEmpty()
    {
        mCobaltData.Ratios.put(chr, List.of(createCobaltRatio(chr, 1000, 0.5, 0.5)));
        mCobaltData.TumorSegments.put(chr, pos(TUMOR_RATIO, 10_000));
        mCobaltData.ReferenceSegments.put(chr, pos(REFERENCE_RATIO, 100_000));

        List<PCFPosition> positions = createPositions();
        assertEquals(2, positions.size());
        checkPcfPosition(positions.get(0), TUMOR_RATIO, 10_000);
        checkPcfPosition(positions.get(1), REFERENCE_RATIO, 100_000);
    }

    @Test
    public void cobaltRatioCommonToTumorAndReference()
    {
        mCobaltData.ReferenceSegments.put(chr, pos(REFERENCE_RATIO, 5000));
        mCobaltData.TumorSegments.put(chr, pos(TUMOR_RATIO, 5000));

        List<PCFPosition> positions = createPositions();
        assertEquals(1, positions.size());

        PCFPosition pcf0 = positions.get(0);
        assertEquals(chr.toString(), pcf0.chromosome());
        assertEquals(5000, pcf0.position());
        //        assertEquals(TUMOR_RATIO, pcf0.Source); // The behaviour is not specified.
    }

    @Test
    public void cobaltEmpty()
    {
        mAmberData.TumorSegments.put(chr, pos(TUMOR_RATIO, 5000));
        mAmberData.TumorSegments.put(chr, pos(TUMOR_RATIO, 6000));
        List<PCFPosition> positions = createPositions();
        assertEquals(2, positions.size());
        checkPcfPosition(positions.get(0), TUMOR_RATIO, 5000);
        checkPcfPosition(positions.get(1), TUMOR_RATIO, 6000);
    }

    @Test
    public void amberAndCobaltHaveData()
    {
        mAmberData.TumorSegments.put(chr, pos(TUMOR_RATIO, 1000));
        mCobaltData.TumorSegments.put(chr, pos(TUMOR_RATIO, 10_000));
        mCobaltData.ReferenceSegments.put(chr, pos(REFERENCE_RATIO, 100_000));

        List<PCFPosition> positions = createPositions();
        assertEquals(3, positions.size());
        checkPcfPosition(positions.get(0), TUMOR_RATIO, 1000);
        checkPcfPosition(positions.get(1), TUMOR_RATIO, 10_000);
        checkPcfPosition(positions.get(2), REFERENCE_RATIO, 100_000);
    }

    void checkPcfPosition(PCFPosition pcfPosition, PCFSource source, int position)
    {
        assertEquals(chr.toString(), pcfPosition.chromosome());
        assertEquals(position, pcfPosition.position());
        assertEquals(source, pcfPosition.Source);
    }

    PCFPosition pos(PCFSource source, int position) {
        return new PCFPosition(source, chr.toString(), position);
    }

    List<PCFPosition> createPositions()
    {
        Map<Chromosome,List<PCFPosition>> map = PCFPositionsSupplier.createPositions(mAmberData, mCobaltData);
        assertEquals(1, map.size());
        return map.get(chr);
    }
}
