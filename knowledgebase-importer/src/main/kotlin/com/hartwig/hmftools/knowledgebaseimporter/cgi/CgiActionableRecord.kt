package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.SequenceVariantType
import com.hartwig.hmftools.knowledgebaseimporter.output.*
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

        operator fun invoke(input: CgiActionableInput, treatmentTypeMap: Map<String, String>): CgiActionableRecord {
            val metadata = CgiMetadata(input.gene, input.transcript ?: "na")
            val events = readSomaticEvents(input)
            val actionability = readActionability(input, treatmentTypeMap)
            if (events.isEmpty()) {
                val aOrBLevelCount = actionability.filter { it.hmfLevel == "A" || it.hmfLevel == "B" }.size
                logger.warn("Could not extract somatic event from:\tcgi\t${input.gene}\t${input.Alteration}\t${input.`Alteration type`}\t$aOrBLevelCount")
            }
            return CgiActionableRecord(metadata, events, actionability, readCgiDrugs(input))
        }

        private fun readSomaticEvents(input: CgiActionableInput): List<SomaticEvent> {
            return listOfNotNull(readProteinAnnotation(input), readCdnaVariant(input), readCNV(input), readFusion(input)) +
                    readGdnaVariants(input) + readGenericMutations(input)
        }

        private fun readAlterations(input: CgiActionableInput) = input.Alteration.substringAfter(":").split(",").map { it.trim() }

        private fun readGdnaVariants(input: CgiActionableInput): List<GDnaVariant> {
            return input.gDNA.orEmpty().split("__").map { it.trim() }.filterNot { it.isBlank() }.map { GDnaVariant(it) }
        }

        private fun readProteinAnnotation(input: CgiActionableInput): ProteinAnnotation? {
            val proteinAnnotation = input.individual_mutation?.substringAfter(':', "")
            return if (proteinAnnotation.isNullOrBlank() || input.transcript.isNullOrBlank()) {
                null
            } else {
                ProteinAnnotation(input.transcript!!, proteinAnnotation!!, SequenceVariantType.OTHER)
            }
        }

        private fun readCdnaVariant(input: CgiActionableInput): CDnaAnnotation? {
            return if (input.cDNA.isNullOrBlank() || input.transcript.isNullOrBlank()) {
                null
            } else {
                CDnaAnnotation(input.transcript!!, input.cDNA!!, SequenceVariantType.OTHER)
            }
        }

        private fun readCNV(input: CgiActionableInput): CnvEvent? = when {
            input.`Alteration type` != "CNA" -> null
            input.Alteration.contains("amp") -> CnvEvent.amplification(input.gene)
            else                             -> CnvEvent.deletion(input.gene)
        }

        private fun readFusion(input: CgiActionableInput): FusionEvent? = when {
            input.`Alteration type` == "FUS"   -> fusionReader.read(input.gene, input.Alteration)
            else                               -> null
        }

        private fun readGenericMutations(input: CgiActionableInput): List<GenericMutation> {
            if (input.`Alteration type` != "MUT") return emptyList()
            return readAlterations(input).mapNotNull {
                when {
                    isAnyMutation(it)   -> GeneMutations(input.gene, input.transcript)
                    isCodonMutation(it) -> CodonMutations(input.gene, input.transcript, codonNumber(it))
                    isCodonRange(it)    -> CodonRangeMutations(input.gene, input.transcript,
                                                               it.substringBefore("-").toInt(), it.substringAfter("-").toInt())
                    else                -> null
                }
            }
        }

        private fun readActionability(input: CgiActionableInput, treatmentTypeMap: Map<String, String>): List<Actionability> {
            val cancerTypes = input.`Primary Tumor type`.split(";").map { it.trim() }
            val level = input.`Evidence level`
            val association = input.Association
            return Actionability("cgi", input.Alteration, cancerTypes, readDrugs(input, treatmentTypeMap), level, association,
                                 "Predictive", HmfLevel(level), HmfResponse(association))
        }

        private fun readDrugs(input: CgiActionableInput, treatmentTypeMap: Map<String, String>): List<HmfDrug> {
            val drugs = readDrugNames(input)
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

        private fun readDrugNames(input: CgiActionableInput): List<String> {
            val drugNames = readDrugsField(input.Drug.orEmpty())
            return if (drugNames.isEmpty()) {
                readDrugsField(input.`Drug family`.orEmpty())
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

        private fun isAnyMutation(alteration: String) = alteration == "."
        private fun isCodonMutation(alteration: String) = alteration.matches(CODON_PATTERN.toRegex(RegexOption.IGNORE_CASE))
        private fun isCodonRange(alteration: String) = alteration.matches("$NUMBER_GROUP_PATTERN-$NUMBER_GROUP_PATTERN".toRegex())


        private fun codonNumber(alteration: String): Int {
            val matchResult = NUMBER_GROUP_PATTERN.toRegex().find(alteration)
            return matchResult!!.groupValues[1].toInt()
        }
    }
}
