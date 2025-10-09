package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltColumns.ENCODED_CHROMOSOME_POS;
import static com.hartwig.hmftools.cobalt.CobaltColumns.POSITION;
import static com.hartwig.hmftools.cobalt.CobaltColumns.RATIO;
import static com.hartwig.hmftools.cobalt.CobaltColumns.READ_GC_CONTENT;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.consolidation.LowCovBucket;
import com.hartwig.hmftools.cobalt.consolidation.LowCoverageRatioMapper;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Before;
import org.junit.Test;

import tech.tablesaw.api.BooleanColumn;
import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.LongColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;

public class BamRatiosTest extends CalculationsTestBase
{

    @Before
    public void setup()
    {
    }

}
