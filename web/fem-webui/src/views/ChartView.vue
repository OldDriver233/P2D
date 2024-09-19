<template>
  <el-row :gutter="15">
    <el-col :span="5">
      <el-card>
        <el-radio-group v-model="selection">
          <el-radio value="u">Concentration of electrolyte</el-radio>
          <el-radio value="delta_u">Voltage delta of 2 ends</el-radio>
          <el-radio value="voltage">Cell voltage</el-radio>
          <el-radio value="c_star">Concentration of particle at surface</el-radio>
        </el-radio-group>
        <el-button @click="update_data">Update</el-button>
      </el-card>
    </el-col>
    <el-col :span="19">
      <el-card>
        <v-chart class="chart" :option="option" />
      </el-card>  
    </el-col>
  </el-row>
</template>
  
<script setup lang="ts">
import { use } from "echarts/core";
import { CanvasRenderer } from "echarts/renderers";
import { LineChart } from "echarts/charts";
import {
  GridComponent
} from "echarts/components";
import VChart, { THEME_KEY } from "vue-echarts";
import { ref, provide, computed, reactive } from "vue";
import axios from "axios";
import type { Ref } from "vue";
import type { ECBasicOption } from "echarts/types/dist/shared";
import { ElRow, ElCol, ElCard, ElMessage, ElRadioGroup, ElRadio } from "element-plus";

use([
  CanvasRenderer,
  LineChart,
  GridComponent
])
  
provide(THEME_KEY, "light")

type ResultType = {
  delta_u: number[],
  voltage: number[],
  c_star: number[],
  u: number[],
}

var result: ResultType = reactive({
  delta_u: [],
  voltage: [],
  c_star: [],
  u: [],
})

const last_update: Ref<number> = ref(0)
const selection: Ref<"u" | "c_star" | "voltage" | "delta_u"> = ref('u')
const x_axis = computed(() => Array.from({length: result[selection.value].length}, (_, i) => i))
const display_data = computed(() => result[selection.value])


function update_data() {
  axios.get("http://localhost:3000/api/result", {headers: {last_update_at: last_update.value}})
    .then(function (resp) {
      last_update.value = resp.data.lastUpdateAt
      let fetched_result: ResultType = resp.data.result
      result.delta_u = fetched_result.delta_u
      result.voltage = fetched_result.voltage
      result.c_star = fetched_result.c_star
      result.u = fetched_result.u
      //result = resp.data.result
      ElMessage({
        message: "Data updated",
        type: "success"
      })
    })
    .catch(function (error) {
      if (error.response) {
        if (error.response.status === 304) {
          ElMessage({
            message: "No new data",
            type: "info"
          })
        } else {
          ElMessage({
            message: "Error: " + error.response.data,
            type: "error"
          })
        }
      } else {
        ElMessage({
          message: "Error: " + error.message,
          type: "error"
        })
      }
    })
}
  
const option: Ref<ECBasicOption> = ref({
  xAxis: {
    data: x_axis
  },
  yAxis: {
    type: 'value'
  },
  series: [
    {
      data: display_data,
      type: 'line',
      smooth: true,
      showSymbol: false
    }
  ],
  itemStyle: {

  }
});
</script>
  
  <style scoped>
  .chart {
    height: 400px;
  }
  </style>