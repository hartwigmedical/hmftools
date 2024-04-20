package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.VariantData;

import htsjdk.variant.variantcontext.VariantContext;

public class RefVariantData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public final VariantType Type;
    public final String Gene;

    public final String CanonicalEffect;
    public final CodingEffect CanonicalCodingEffect;
    public final String HgvsCodingImpact;
    public final String HgvsProteinImpact;

    public final CodingEffect WorstCodingEffect;

    public final int LocalPhaseSet;
    public final String Microhomology;
    public final String RepeatSequence;
    public final int RepeatCount;
    public final boolean Reported;
    public final boolean IsHotspot;

    private Integer mLiftedPosition; // if liftover has been run

    public RefVariantData(
            final String chromosome, final int position, final String ref, final String alt, final VariantType type,
            final String gene, final String canonicalEffect, final CodingEffect canonicalCodingEffect,
            final CodingEffect worstCodingEffect, final String hgvsCodingImpact, final String hgvsProteinImpact,
            final String microhomology, final String repeatSequence, int repeatCount, int localPhaseSet,
            boolean reported, boolean isHotspot)
    {
        Chromosome = chromosome;
        Position = position;
        mLiftedPosition = null;
        Ref = ref;
        Alt = alt;
        Type = type;
        Gene = gene;

        CanonicalEffect = canonicalEffect;
        CanonicalCodingEffect = canonicalCodingEffect;
        WorstCodingEffect = worstCodingEffect;

        LocalPhaseSet = localPhaseSet;
        HgvsCodingImpact = hgvsCodingImpact;
        HgvsProteinImpact = hgvsProteinImpact;
        Microhomology = microhomology;
        RepeatSequence = repeatSequence;
        RepeatCount = repeatCount;
        Reported = reported;
        IsHotspot = isHotspot;
    }

    public void setLiftedPosition(final int position) { mLiftedPosition = position; }
    public Integer liftedPosition() { return mLiftedPosition; }

    public static RefVariantData fromSomatic(final SomaticVariant variant)
    {
        int localPhaseSet = variant.topLocalPhaseSet() != null ? variant.topLocalPhaseSet() : NO_LOCAL_PHASE_SET;

        return new RefVariantData(
                variant.chromosome(), (int)variant.position(), variant.ref(), variant.alt(), variant.type(), variant.gene(),
                variant.canonicalEffect(), variant.canonicalCodingEffect(), variant.worstCodingEffect(),
                variant.canonicalHgvsCodingImpact(), variant.canonicalHgvsProteinImpact(),
                variant.microhomology(), variant.repeatSequence(), variant.repeatCount(),
                localPhaseSet, variant.reported(), variant.isHotspot());
    }

    public static List<RefVariantData> loadVariantsFromVcf(final String filename, final ComparisonConfig config)
    {
        VcfFileReader vcfFileReader = new VcfFileReader(filename);

        if(!vcfFileReader.fileValid())
        {
            PV_LOGGER.error("failed to read VCF({})", filename);
            return null;
        }

        List<RefVariantData> variants = Lists.newArrayList();

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(variantContext.isFiltered())
                continue;
            VariantContextDecorator variant = new VariantContextDecorator(variantContext);

            String chromosome = variant.chromosome();
            int position = variant.position();
            Integer liftedPosition = null;

            if(config.LiftoverCache != null)
            {
                liftedPosition = position;
                position = config.LiftoverCache.convertPositionTo38(chromosome, position);
                chromosome = config.RefGenVersion.versionedChromosome(chromosome);
            }

            VariantImpact impact = variant.variantImpact();

            int localPhaseSet = variant.localPhaseSet() != null ? variant.localPhaseSet() : NO_LOCAL_PHASE_SET;

            RefVariantData refVariant = new RefVariantData(
                    chromosome, position, variant.ref(), variant.alt(), variant.type(), variant.gene(),
                    impact.CanonicalEffect, impact.CanonicalCodingEffect, impact.WorstCodingEffect,
                    impact.CanonicalHgvsCoding, impact.CanonicalHgvsProtein,
                    variant.microhomology(), variant.repeatSequence(), variant.repeatCount(),
                    localPhaseSet, variant.reported(), variant.isHotspot());

            if(liftedPosition != null)
                refVariant.setLiftedPosition(liftedPosition);

            variants.add(refVariant);
        }

        return variants;
    }

    public boolean matches(final VariantData variant)
    {
        return variant.Chromosome.equals(Chromosome) && variant.Position == Position && variant.Ref.equals(Ref) && variant.Alt.equals(Alt);
    }

    public String toString()
    {
        return String.format("pos(%s:%d) variant(%s: %s>%s) canon(%s: %s) worst(%s) hgvs(coding=%s protein=%s)",
                Chromosome, Position, Type, Ref, Alt, CanonicalCodingEffect, CanonicalEffect,
                WorstCodingEffect, HgvsCodingImpact, HgvsProteinImpact);
    }

    public static String tsvHeader()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("sampleId");
        sj.add("chromosome");
        sj.add("position");
        sj.add("type");
        sj.add("ref");
        sj.add("alt");
        sj.add("gene");
        sj.add("worstCodingEffect");
        sj.add("canonicalEffect");
        sj.add("canonicalCodingEffect");
        sj.add("canonicalHgvsCodingImpact");
        sj.add("canonicalHgvsProteinImpact");
        sj.add("microhomology");
        sj.add("repeatSequence");
        sj.add("repeatCount");
        sj.add("localPhaseSet");
        sj.add("reported");
        sj.add("hotspot");
        return sj.toString();
    }

    public String tsvData(final String sampleId)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(sampleId);
        sj.add(Chromosome);
        sj.add(String.valueOf(Position));
        sj.add(Type.toString());
        sj.add(Ref);
        sj.add(Alt);
        sj.add(Gene);
        sj.add(String.valueOf(WorstCodingEffect));
        sj.add(CanonicalEffect);
        sj.add(String.valueOf(CanonicalCodingEffect));
        sj.add(HgvsCodingImpact);
        sj.add(HgvsProteinImpact);
        sj.add(Microhomology);
        sj.add(RepeatSequence);
        sj.add(String.valueOf(RepeatCount));
        sj.add(LocalPhaseSet == NO_LOCAL_PHASE_SET ? "NULL" : String.valueOf(LocalPhaseSet));
        sj.add(Reported ? "1" : "0");
        sj.add(IsHotspot ? Hotspot.HOTSPOT.toString() : "NONE");
        return sj.toString();
    }
}
