
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - Gene Reference Data Utilities"

dependencies {
    implementation(project(":hmf-common"))
    implementation(libs.jooq)
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.geneutils.CoOccurenceCalcs")
}
