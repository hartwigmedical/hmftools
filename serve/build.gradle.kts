
plugins {
    id("com.hartwig.java-conventions")
    id("com.hartwig.tests-jar") // create a jar for tests
    application
}

description = "HMF Tools - SERVE"

dependencies {
    api(project(":hmf-common"))
    implementation(project(":ckb-importer"))
    implementation(project(":iclusion-importer"))
    implementation(project(":vicc-importer"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.serve.ServeApplication")
}
