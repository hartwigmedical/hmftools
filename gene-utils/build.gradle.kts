
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - Gene Reference Data Utilities"

dependencies {
    implementation(project(":hmf-common"))
    implementation(libs.jooq)
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.geneutils.CoOccurenceCalcs")
}
