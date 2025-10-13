package com.hartwig.hmftools.cobalt.diploid;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;

public class DiploidRegionLoader implements Consumer<Locatable>
{
    private final Table mResult;
    private final Table mContigResult;
    private final ListMultimap<Chromosome, DiploidStatus> Regions = ArrayListMultimap.create();


    private final ChromosomePositionCodec mChromosomePosCodec;

    private String mChromosome = null;
    private Chromosome CurrentChromosome;
    private int mPosition = 0;

    public DiploidRegionLoader(final ChromosomePositionCodec chromosomePosCodec, final String diploidBedPath) throws IOException
    {
        mChromosomePosCodec = chromosomePosCodec;
        mResult = Table.create(StringColumn.create(CobaltColumns.CHROMOSOME), IntColumn.create(CobaltColumns.POSITION));
        mContigResult = mResult.emptyCopy();

        if(diploidBedPath != null)
        {
            List<BEDFeature> bedFeatures = new ArrayList<>();

            CB_LOGGER.info("Reading diploid regions from {}", diploidBedPath);
            try(final AbstractFeatureReader<BEDFeature, LineIterator> reader = getFeatureReader(diploidBedPath,
                    new BEDCodec(),
                    false))
            {
                for(BEDFeature bedFeature : reader.iterator())
                {
                    bedFeatures.add(bedFeature);
                }
            }

            bedFeatures.forEach(this);
        }
    }

    @Override
    public void accept(@NotNull Locatable bed)
    {
        if(mChromosome == null || !bed.getContig().equals(mChromosome))
        {
            finaliseCurrent();
            mChromosome = bed.getContig();
            CurrentChromosome = HumanChromosome.fromString(mChromosome);
            mPosition = 1;
        }
        createRatio(bed.getContig(), bed.getStart(), bed.getEnd());
        mPosition = bed.getEnd() + 1;
    }

    private void createRatio(String contig, int start, int end)
    {
        // Add non-diploid entries from current position to the start of this block.
        int nonDiploidPosition = mPosition;
        while(nonDiploidPosition < start)
        {
            DiploidStatus region = new DiploidStatus(contig, nonDiploidPosition, nonDiploidPosition + WINDOW_SIZE - 1, false);
            Regions.put(CurrentChromosome, region);
            nonDiploidPosition += WINDOW_SIZE;
        }

        int position = start;
        while(position < end)
        {
            DiploidStatus region = new DiploidStatus(contig, position, position + WINDOW_SIZE - 1, true);
            Regions.put(CurrentChromosome, region);

            Row row = mContigResult.appendRow();
            row.setString(CobaltColumns.CHROMOSOME, contig);
            row.setInt(CobaltColumns.POSITION, position);
            position += WINDOW_SIZE;
        }
    }

    private void finaliseCurrent()
    {
        if(mChromosome != null && mPosition > 0)
        {
            mResult.append(mContigResult);
        }

        mContigResult.clear();
    }

    @NotNull
    public Table build()
    {
        finaliseCurrent();
        mChromosomePosCodec.addEncodedChrPosColumn(mResult, false);
        return mResult;
    }

    public ListMultimap<Chromosome, DiploidStatus> regions()
    {
        return Regions;
    }
}
