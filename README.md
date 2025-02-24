# Beginner guide and toy code implementation of Fully Homomorphic Encryption : FHEW (Korean)

동형암호는 암호화 된 상태에서 연산을 수행할 수 있는 암호 시스템입니다. 지원하는 연산에 따라 다양한 동형암호들이 존재하는데 크게 다음과 같습니다.

1. 실수에 대한 연산을 제공하는 <a href = "https://eprint.iacr.org/2016/421">CKKS</a>.
2. 정수에 대한 연산을 제공하는 <a href = "https://eprint.iacr.org/2011/277">BGV</a> 와 <a href = "https://eprint.iacr.org/2012/144">BFV</a>.
3. 비트에 대한 논리연산을 제공하는 <a href = "https://eprint.iacr.org/2014/816">FHEW</a>/<a href = "https://eprint.iacr.org/2016/870">TFHE</a>.

여기서는 FHEW 에 대한 내용을 다루고 있으며, 실제로 Toy code 까지 만들어나가게 됩니다.

참고할 논문은 다음 두 가지 입니다.

* <a href="https://eprint.iacr.org/2014/816">FHEW: Bootstrapping Homomorphic Encryption in less than a second</a>
* <a href="https://eprint.iacr.org/2020/086">Bootstrapping in FHEW-like Cryptosystems</a>

크게 다음과 같은 내용들을 다루고 코드로 구현해봅니다.

1. 대수적 구조체 $\mathbb{Z}_q[X]/(X^N+1)$

2. Learning with errors problem (LWE) 와 LWE 에서 정의되는 연산들.

3. Ring Learning with errors problem (RLWE) 와 RLWE 에서 정의되는 연산들.

4. RLWE 의 변형인 RLWE' 과 RGSW.

5. Key switching, Modulus Switching.

6. Blind rotation.

---

> 혹시 질문이 있으시거나 내용에 오류가 있으면 말씀해주시기 바랍니다.
> deokhwa@inha.edu
