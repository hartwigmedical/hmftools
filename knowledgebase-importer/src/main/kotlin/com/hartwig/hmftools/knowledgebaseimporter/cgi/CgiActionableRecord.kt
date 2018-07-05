package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager

data class CgiActionableRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                               override val actionability: List<Actionability>, val cgiDrugs: List<CgiDrug>) : RecordMetadata by metadata,
        ActionableRecord {
    companion object {
        private val logger = LogManager.getLogger("CgiActionableRecord")
        private const val NUMBER_GROUP_PATTERN = "([0-9]+)"
        private const val AMINO_ACID_LETTERS = "GALMFWKQESPVICYHRNDT"
        private const val AMINO_ACID_CODON_PATTERN = "[$AMINO_ACID_LETTERS]$NUMBER_GROUP_PATTERN\\."
        private const val ANY_CODON_PATTERN = "\\.$NUMBER_GROUP_PATTERN\\."
        private const val CODON_PATTERN = "$AMINO_ACID_CODON_PATTERN|$ANY_CODON_PATTERN"
        private val FUSION_SEPARATORS = setOf("__")
        private val FUSIONS_TO_FLIP = setOf(FusionPair("ABL1", "BCR"),
                                            FusionPair("PDGFRA", "FIP1L1"),
                                            FusionPair("PDGFB", "COL1A1"))
        private val FUSIONS_TO_FILTER = setOf(FusionPair("RET", "TPCN1"))
        private val fusionReader = FusionReader(separators = FUSION_SEPARATORS, filterSet = FUSIONS_TO_FILTER, flipSet = FUSIONS_TO_FLIP)

        operator fun invoke(record: CSVRecord, treatmentTypeMap: Map<String, String>): CgiActionableRecord {
            val metadata = CgiMetadata(record["Gene"], record["transcript"] ?: "na")
            val events = readSomaticEvents(record)
            if (events.isEmpty()) {
                logger.warn("Could not extract somatic event from:\tcgi\t${record["Gene"]}\t${record["Alteration"]}\t${record["Alteration type"]}")
            }
            val actionability = readActionability(record, treatmentTypeMap).filterNot { it.significance == "No Responsive" }
            return CgiActionableRecord(metadata, events, actionability, readCgiDrugs(record))
        }

        private fun readSomaticEvents(record: CSVRecord): List<SomaticEvent> {
            return listOfNotNull(readProteinAnnotation(record), readCdnaVariant(record), readCNV(record), readFusion(record)) +
                    readGdnaVariants(record) + readGenericMutations(record)
        }

        private fun readAlterations(record: CSVRecord) = record["Alteration"].substringAfter(":").split(",").map { it.trim() }

        private fun readGdnaVariants(record: CSVRecord): List<GDnaVariant> {
            return record["gDNA"].orEmpty().split("__").map { it.trim() }.filterNot { it.isBlank() }.map { GDnaVariant(it) }
        }

        private fun readProteinAnnotation(record: CSVRecord): ProteinAnnotation? {
            val transcript = record["transcript"]
            val proteinAnnotation = record["individual_mutation"]?.substringAfter(':', "")
            return if (proteinAnnotation.isNullOrBlank() || transcript.isNullOrBlank()) {
                null
            } else {
                ProteinAnnotation(transcript, proteinAnnotation!!)
            }
        }

        private fun readCdnaVariant(record: CSVRecord): CDnaAnnotation? {
            val transcript = record["transcript"]
            val cdnaAnnotation = record["cDNA"]
            return if (cdnaAnnotation.isNullOrBlank() || transcript.isNullOrBlank()) {
                null
            } else {
                CDnaAnnotation(transcript, cdnaAnnotation!!)
            }
        }

        private fun readCNV(record: CSVRecord): CnvEvent? = when {
            record["Alteration type"] != "CNA"   -> null
            record["Alteration"].contains("amp") -> CnvEvent(record["Gene"], "Amplification")
            else                                 -> CnvEvent(record["Gene"], "Deletion")
        }

        private fun readFusion(record: CSVRecord): FusionEvent? = when {
            record["Alteration type"] == "FUS" -> fusionReader.read(record["Gene"], record["Alteration"])
            else                               -> null
        }

        private fun readGenericMutations(record: CSVRecord): List<GenericMutation> {
            if (record["Alteration type"] != "MUT") return emptyList()
            return readAlterations(record).mapNotNull {
                when {
                    isGeneMutation(it)  -> GeneMutations(record["Gene"], record["transcript"])
                    isCodonMutation(it) -> CodonMutations(record["Gene"], record["transcript"], codonNumber(it))
                    isCodonRange(it)    -> CodonRangeMutations(record["Gene"], record["transcript"],
                                                               it.substringBefore("-").toInt(), it.substringAfter("-").toInt())
                    else                -> null
                }
            }
        }

        private fun readActionability(record: CSVRecord, treatmentTypeMap: Map<String, String>): List<Actionability> {
            val cancerTypes = record["Primary Tumor type"].split(";").map { it.trim() }
            val level = record["Evidence level"]
            val association = record["Association"]
            return Actionability("cgi", record["Alteration"], cancerTypes, readDrugs(record, treatmentTypeMap), level, association,
                                 "Predictive", HmfLevel(level), HmfResponse(association))
        }

        private fun readDrugs(record: CSVRecord, treatmentTypeMap: Map<String, String>): List<HmfDrug> {
            val drugs = readDrugNames(record)
            return drugs.map { name ->
                if (name.contains("+")) {
                    val type = name.split("+").map { it.trim() }
                            .joinToString(" + ") { treatmentTypeMap[it.toLowerCase()] ?: "Unknown" }
                    HmfDrug(name, type)
                } else {
                    HmfDrug(name, treatmentTypeMap[name.toLowerCase()] ?: "Unknown")
                }
            }
        }

        private fun readDrugNames(record: CSVRecord): List<String> {
            val drugNames = readDrugsField(record["Drug"].orEmpty())
            return if (drugNames.isEmpty()) {
                readDrugsField(record["Drug family"].orEmpty())
            } else {
                drugNames
            }
        }

        private fun readDrugsField(drugField: String): List<String> {
            return drugField.replace(";", " + ")
                    .removeSurrounding("[", "]")
                    .split(",")
                    .filterNot { it.isBlank() }
        }

        private fun isGeneMutation(alteration: String) = alteration == "."
        private fun isCodonMutation(alteration: String) = alteration.matches(CODON_PATTERN.toRegex(RegexOption.IGNORE_CASE))
        private fun isCodonRange(alteration: String) = alteration.matches("$NUMBER_GROUP_PATTERN-$NUMBER_GROUP_PATTERN".toRegex())


        private fun codonNumber(alteration: String): Int {
            val matchResult = NUMBER_GROUP_PATTERN.toRegex().find(alteration)
            return matchResult!!.groupValues[1].toInt()
        }
    }
}
