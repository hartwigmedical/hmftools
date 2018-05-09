package com.hartwig.hmftools.knowledgebaseimporter.output

interface ActionableItem<out T> {
    val actionability: Actionability
    val event: T
}