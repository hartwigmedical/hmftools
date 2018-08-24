package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

import org.apache.logging.log4j.LogManager

interface CorrectedInput<T : CorrectedInput<T>> {
    companion object {
        private val logger = LogManager.getLogger(this::class.java)
    }

    fun correct(): T

    fun corrected(): T {
        val correctedValue = correct()
        if (this != correctedValue) logger.info("Corrected $this to $correctedValue")
        return correctedValue
    }
}
