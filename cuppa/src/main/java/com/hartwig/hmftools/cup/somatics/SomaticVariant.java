package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_FLAG;
import static com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser.VAR_IMPACT;
import static com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser.fromVariantContext;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;

import org.jooq.Record;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariant
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public final VariantType Type;
    public final String Gene;
    public final String TrinucleotideContext;
    public final int RepeatCount;

    public static final String FLD_CHR = "Chromosome";
    public static final String FLD_POSITION = "Position";
    public static final String FLD_REF = "Ref";
    public static final String FLD_ALT = "Alt";
    public static final String FLD_TYPE = "Type";
    public static final String FLD_REPEAT_COUNT = "RepeatCount";
    public static final String FLD_GENE = "Gene";
    public static final String FLD_TRINUC_CONTEXT = "TriNucContext";

    public static final String SOMATIC_VAR_FILE_ID = ".somatic_variants.csv";
    public static final String SOMATIC_VARIANTS_DIR_CFG = "somatic_variants_dir";
    public static final String SOMATIC_VARIANTS_DIR_DESC = "Directory containing somatic variant generic tabular files";

    public SomaticVariant(
            final String chromosome, final int position, final String ref, final String alt, final VariantType type,
            final String gene, final String trinucleotideContext, final int repeatCount)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        Gene = gene;
        TrinucleotideContext = trinucleotideContext;
        RepeatCount = repeatCount;
    }

    public static String generateFilename(final String baseDir, final String sampleId)
    {
        return FileWriterUtils.checkAddDirSeparator(baseDir) + sampleId + SOMATIC_VAR_FILE_ID;
    }

    public static SomaticVariant fromContext(final VariantContext variantContext)
    {
        int position = variantContext.getStart();
        String chromosome = variantContext.getContig();

        String ref = variantContext.getReference().getBaseString();
        String alt = !variantContext.getAlternateAlleles().isEmpty() ? variantContext.getAlternateAlleles().get(0).toString() : ref;

        if(alt.equals("*") || alt.equals("N")) // unhandled for now
            alt = ref;

        String gene = "";

        if(variantContext.hasAttribute(VAR_IMPACT))
        {
            VariantImpact impact = fromVariantContext(variantContext);
            gene = impact.CanonicalGeneName;
        }

        return new SomaticVariant(
                chromosome, position, ref, alt, VariantType.type(variantContext), gene,
                variantContext.getAttributeAsString(TRINUCLEOTIDE_FLAG, ""),
                variantContext.getAttributeAsInt(REPEAT_COUNT_FLAG, 0));
    }

    public static SomaticVariant fromRecord(final Record record)
    {
        return new SomaticVariant(
                record.getValue(Tables.SOMATICVARIANT.CHROMOSOME),
                record.getValue(Tables.SOMATICVARIANT.POSITION),
                record.getValue(Tables.SOMATICVARIANT.REF),
                record.getValue(Tables.SOMATICVARIANT.ALT),
                VariantType.valueOf(record.getValue(SOMATICVARIANT.TYPE)),
                record.getValue(Tables.SOMATICVARIANT.GENE),
                record.getValue(SOMATICVARIANT.TRINUCLEOTIDECONTEXT),
                0); // features makes its own DB call for specific INDELs so isn't retrieved here
    }

    public static String csvHeader()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(FLD_CHR);
        sj.add(FLD_POSITION);
        sj.add(FLD_REF);
        sj.add(FLD_ALT);
        sj.add(FLD_TYPE);
        sj.add(FLD_GENE);
        sj.add(FLD_TRINUC_CONTEXT);
        sj.add(FLD_REPEAT_COUNT);
        return sj.toString();
    }

    public String toString()
    {
        return String.format("pos(%s:%d) variant(%s: %s>%s) gene(%s) triNuc(%s)",
                Chromosome, Position, Type, Ref, Alt, Gene, TrinucleotideContext);
    }
}
