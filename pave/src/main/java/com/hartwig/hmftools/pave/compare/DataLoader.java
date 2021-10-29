package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.pave.GeneDataCache;

import org.jooq.Record;
import org.jooq.Record20;
import org.jooq.Result;

public final class DataLoader
{
    public static List<RefVariantData> loadSampleDatabaseRecords(
            final String sampleId, final DatabaseAccess dbAccess, final GeneDataCache geneDataCache)
    {
        List<RefVariantData> variants = Lists.newArrayList();

        Result<Record20<String, String, Integer, String, String, String, String, String, String, String, Integer, String, String, String,
                String, Integer, Integer, Integer, Byte, String>>
                result = dbAccess.context()
                .select(SOMATICVARIANT.GENE, SOMATICVARIANT.CHROMOSOME, SOMATICVARIANT.POSITION,
                        SOMATICVARIANT.REF, SOMATICVARIANT.ALT, SOMATICVARIANT.TYPE, SOMATICVARIANT.GENE,
                        SOMATICVARIANT.CANONICALEFFECT, SOMATICVARIANT.CANONICALCODINGEFFECT, SOMATICVARIANT.WORSTCODINGEFFECT,
                        SOMATICVARIANT.GENESAFFECTED, SOMATICVARIANT.CANONICALHGVSCODINGIMPACT, SOMATICVARIANT.CANONICALHGVSPROTEINIMPACT,
                        SOMATICVARIANT.MICROHOMOLOGY, SOMATICVARIANT.REPEATSEQUENCE, SOMATICVARIANT.REPEATCOUNT,
                        SOMATICVARIANT.PHASEDINFRAMEINDEL, SOMATICVARIANT.LOCALPHASESET, SOMATICVARIANT.REPORTED, SOMATICVARIANT.HOTSPOT)
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.GENE.isNotNull())
                .orderBy(SOMATICVARIANT.CHROMOSOME, SOMATICVARIANT.POSITION)
                .fetch();

        for(Record record : result)
        {
            //final SomaticVariant variant = SomaticVariantDAO.buildFromRecord(record);

            if(geneDataCache != null)
            {
                String gene = record.getValue(SOMATICVARIANT.GENE);

                if(!geneDataCache.isDriverPanelGene(gene))
                    continue;
            }

            variants.add(RefVariantData.fromRecord(record));
        }

        return variants;
    }

    public static Map<String,List<RefVariantData>> processRefVariantFile(final String filename)
    {
        Map<String,List<RefVariantData>> sampleVariantsMap = Maps.newHashMap();

        if(filename == null)
            return sampleVariantsMap;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            final String fileDelim = "\t";
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, fileDelim);

            int sampleIndex = fieldsIndexMap.get("sampleId");

            int chrIndex = fieldsIndexMap.get("chromosome");
            int posIndex = fieldsIndexMap.get("position");
            int refIndex = fieldsIndexMap.get("ref");
            int altIndex = fieldsIndexMap.get("alt");
            int typeIndex = fieldsIndexMap.get("type");
            int geneIndex = fieldsIndexMap.get("gene");
            int canonicalEffectIndex = fieldsIndexMap.get("canonicalEffect");
            int canonicalCodingEffectIndex = fieldsIndexMap.get("canonicalCodingEffect");
            int worstCodingEffectIndex = fieldsIndexMap.get("worstCodingEffect");
            int genesAffectedIndex = fieldsIndexMap.get("genesEffected");
            int canonicalHgvsCodingImpactIndex = fieldsIndexMap.get("canonicalHgvsCodingImpact");
            int canonicalHgvsProteinImpactIndex = fieldsIndexMap.get("canonicalHgvsProteinImpact");
            int microhomologyIndex = fieldsIndexMap.get("microhomology");
            int repeatSequenceIndex = fieldsIndexMap.get("repeatSequence");
            int repeatCountIndex = fieldsIndexMap.get("repeatCount");
            int phasedInframeIndelIndex = fieldsIndexMap.get("phasedInframeIndel");
            int localPhaseSetIndex = fieldsIndexMap.get("localPhaseSet");
            int reportedIndex = fieldsIndexMap.get("reported");
            Integer hotspotIndex = fieldsIndexMap.get("hotspot");

            List<RefVariantData> variants = null;

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(fileDelim, -1);

                String sampleId = items[sampleIndex];

                variants = sampleVariantsMap.get(sampleId);
                if(variants == null)
                {
                    variants = Lists.newArrayList();
                    sampleVariantsMap.put(sampleId, variants);
                }

                int localPhaseSet = !items[localPhaseSetIndex].equals("NULL") ? Integer.parseInt(items[localPhaseSetIndex]) : NO_LOCAL_PHASE_SET;

                RefVariantData variant = new RefVariantData(
                        items[chrIndex], Integer.parseInt(items[posIndex]), items[refIndex], items[altIndex],
                        VariantType.valueOf(items[typeIndex]), items[geneIndex], items[canonicalEffectIndex],
                        items[canonicalCodingEffectIndex].isEmpty() ? NONE : CodingEffect.valueOf(items[canonicalCodingEffectIndex]),
                        CodingEffect.valueOf(items[worstCodingEffectIndex]), Integer.parseInt(items[genesAffectedIndex]),
                        items[canonicalHgvsCodingImpactIndex], items[canonicalHgvsProteinImpactIndex],
                        items[microhomologyIndex], items[repeatSequenceIndex], Integer.parseInt(items[repeatCountIndex]),
                        Boolean.parseBoolean(items[phasedInframeIndelIndex]),
                        localPhaseSet, items[reportedIndex].equals("1"),
                        hotspotIndex != null ? items[hotspotIndex].equals(Hotspot.HOTSPOT.toString()) : false);

                variants.add(variant);
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to read ref variant data file: {}", e.toString());
        }

        return sampleVariantsMap;
    }

}
