package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.FLD_ALT;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.FLD_CHR;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.FLD_GENE;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.FLD_POSITION;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.FLD_REF;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.FLD_REPEAT_COUNT;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.FLD_TRINUC_CONTEXT;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.FLD_TYPE;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.file.FileReaderUtils;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.cup.prep.PrepConfig;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class SomaticVariantsLoader
{
    public static final String SOMATIC_VARIANTS_DIR = "somatic_variants_dir";

    public static List<SomaticVariant> loadFromConfig(
            final PrepConfig config, final String sampleId, @Nullable final List<VariantType> variantTypes) throws NoSuchFileException
    {
        File genericVariantsFile = new File(config.somaticVariantsGenericFile(sampleId));
        File vcfFile = new File(config.purpleSomaticVcfFile(sampleId));

        List<SomaticVariant> variants;
        if(genericVariantsFile.isFile())
        {
            if(vcfFile.isFile())
            {
                CUP_LOGGER.error("VCF and generic variants files both exist for sample({})", sampleId);
            }

            variants = loadFromGenericFile(genericVariantsFile.getAbsolutePath(), variantTypes);
        }
        else if(vcfFile.isFile())
        {
            variants = loadFromVcf(vcfFile.getAbsolutePath(), variantTypes);
        }
        else
        {
            throw new NoSuchFileException(String.format("%s or %s not provided", PURPLE_DIR_CFG, SOMATIC_VARIANTS_DIR));
        }

        return variants;
    }

    private static List<SomaticVariant> loadFromVcf(final String vcfFile, @Nullable final List<VariantType> variantTypes)
    {
        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        List<SomaticVariant> variants = new ArrayList<>();

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        if(!vcfFileReader.fileValid())
        {
            CUP_LOGGER.error("Invalid somatic VCF file({})", vcfFile);
            System.exit(1);
        }

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(!filter.test(variantContext))
                continue;

            if(variantTypes == null || variantTypes.contains(VariantType.type(variantContext)))
            {
                variants.add(SomaticVariant.fromContext(variantContext));
            }
        }

        return variants;
    }

    private static List<SomaticVariant> loadFromGenericFile(final String filename, @Nullable final List<VariantType> variantTypes)
    {
        List<SomaticVariant> variants = new ArrayList<>();

        if(filename == null || filename.isEmpty())
            System.exit(1);

        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = FileReaderUtils.createFieldsIndexMap(header, CSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHR);
            int posIndex = fieldsIndexMap.get(FLD_POSITION);
            int refIndex = fieldsIndexMap.get(FLD_REF);
            int altIndex = fieldsIndexMap.get(FLD_ALT);
            int typeIndex = fieldsIndexMap.get(FLD_TYPE);
            int rcIndex = fieldsIndexMap.get(FLD_REPEAT_COUNT);
            int tnIndex = fieldsIndexMap.get(FLD_TRINUC_CONTEXT);
            int geneIndex = fieldsIndexMap.get(FLD_GENE);

            for(final String line : lines)
            {
                final String[] values = line.split(CSV_DELIM);
                SomaticVariant variant = new SomaticVariant(
                        values[chrIndex], Integer.parseInt(values[posIndex]), values[refIndex], values[altIndex],
                        VariantType.valueOf(values[typeIndex]), values[geneIndex], values[tnIndex], Integer.parseInt(values[rcIndex]));

                if(variantTypes == null && !variantTypes.contains(variant.Type))
                    continue;

                variants.add(variant);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read somatic variant flat file({}): {}", filename, e.toString());
            System.exit(1);
        }

        return variants;
    }
}
