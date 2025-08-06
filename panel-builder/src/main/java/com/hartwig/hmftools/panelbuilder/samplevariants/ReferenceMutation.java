package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildMutationProbe;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.wisp.CategoryType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ReferenceMutation extends Variant
{
    private final BasePosition mPosition;
    private final String mRef;
    private final String mAlt;
    private final String mSource;

    private static final Logger LOGGER = LogManager.getLogger(ReferenceMutation.class);

    public ReferenceMutation(
            final BasePosition position, final String ref, final String alt, final String source)
    {
        mPosition = position;
        mRef = ref;
        mAlt = alt;
        mSource = source;
    }

    @Override
    public CategoryType categoryType()
    {
        return CategoryType.REFERENCE;
    }

    @Override
    public String description()
    {
        if(!mRef.isEmpty() && !mAlt.isEmpty())
        {
            return format("%s %s>%s %s", mPosition, mRef, mAlt, mSource);
        }
        else
        {
            return format("%s %s", mPosition, mSource);
        }
    }

    @Override
    public double copyNumber()
    {
        return 0;
    }

    @Override
    public double vaf()
    {
        return 0;
    }

    @Override
    public int tumorFragments()
    {
        return 0;
    }

    @Override
    public boolean reported()
    {
        return false;
    }

    @Override
    public VariantProbeData generateProbe(final RefGenomeInterface refGenome)
    {
        if(mAlt.isEmpty())
        {
            ChrBaseRegion region = probeRegionCenteredAt(mPosition.Chromosome, mPosition.Position);
            return new VariantProbeData(null, region, null, null);
        }
        else
        {
            return buildMutationProbe(mPosition.Chromosome, mPosition.Position, mRef, mAlt, PROBE_LENGTH, refGenome);
        }
    }

    @Override
    public boolean checkFilters()
    {
        return false;
    }

    @Override
    public boolean checkAndRegisterLocation(ProximateLocations registeredLocations)
    {
        registeredLocations.addRegisteredLocation(mPosition.Chromosome, mPosition.Position);
        return true;
    }

    public String toString()
    {
        return description();
    }

    public static List<Variant> loadKnownMutations(final String filename)
    {
        List<Variant> variants = new ArrayList<>();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);

            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posIndex = fieldsIndexMap.get(FLD_POSITION);
            int refIndex = fieldsIndexMap.get(FLD_REF);
            int altIndex = fieldsIndexMap.get(FLD_ALT);
            int sourceIndex = fieldsIndexMap.get("Source");

            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                variants.add(new ReferenceMutation(
                        new BasePosition(values[chrIndex], Integer.parseInt(values[posIndex])),
                        values[refIndex], values[altIndex],
                        values[sourceIndex]));
            }

            LOGGER.info("loaded {} reference variants from file({})", variants.size(), filename);
        }
        catch(Exception e)
        {
            LOGGER.error("failed to load reference variants file({}): {}", filename, e.toString());
            return null;
        }

        return variants;
    }
}
