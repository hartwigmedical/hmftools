
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - Cuppa"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    implementation(libs.jooq)
    implementation(libs.commons.lang3)
    testImplementation(libs.junit)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
}

description = "HMF Tools - Cuppa"

application {
    mainClass.set("com.hartwig.hmftools.cup.CupAnalyser")
}
