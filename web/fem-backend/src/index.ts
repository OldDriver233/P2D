import fastify from 'fastify'
import Fastify from 'fastify'
import {createClient} from 'redis'


const startServer = async () => {
    const client = createClient()
    client.on('error', (err) => console.log('Redis Client Error', err))

    await client.connect()

    const server = Fastify({
        logger: true
    })

    server.get('/api/result', async (request, reply) => {
        const raw_result = await client.lRange('us', 0, -1)
        const result = raw_result.map((item) => Number.parseFloat(item))
        reply.header('Access-Control-Allow-Origin', '*')
        reply.send(result)
    })

    server.listen({port: 3000}, function (err, _) {
        if (err) {
            server.log.error(err)
            process.exit(1)
        }
    })
}

startServer()
